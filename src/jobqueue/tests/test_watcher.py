"""The ESO watcher against a fake TAP: stability, the transfer-log match, partial nights, no-data notes, Artemis."""
import datetime as dt
import os

from jobqueue.tests.fakes import PROGS, FakeEso, write_fits
from jobqueue.watcher import Watcher

NIGHT = '20261008'


def at(h, m=0, day=9):
    return dt.datetime(2026, 10, day, h, m)


def poll(cfg, store, eso, now):
    return Watcher(cfg, store, eso, now=now, log=lambda m: None).run()


def p0_jobs(store, tel=None):
    return [j for j in store.jobs(classes=['P0']) if tel is None or j['telescope'] == tel]


def transfer_log(cfg, tel, lines):
    d = os.path.join(cfg['basedir'], 'Observations', tel)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, 'transfer_log.txt'), 'w') as f:
        f.write(''.join('{} {}\n'.format(n, c) for n, c in lines))


def only(cfg, *tels):
    return dict(cfg, telescopes={t: cfg['telescopes'][t] for t in tels})


def test_watches_last_night_and_the_one_before(cfg, store):
    eso = FakeEso()
    poll(only(cfg, 'Europa'), store, eso, at(14))
    assert [c[1] for c in eso.calls if c[0] == 'summary'] == ['20261008', '20261007']
    eso.calls.clear()
    poll(only(cfg, 'Europa'), store, eso, at(10, 30))       # before 11:00 UTC the last night is still the 7th
    assert [c[1] for c in eso.calls if c[0] == 'summary'] == ['20261007', '20261006']


def test_stable_count_queues_download_then_pipeline_once(cfg, store):
    cfg = only(cfg, 'Europa')
    eso = FakeEso()
    eso.objects[(PROGS['Europa'], NIGHT)] = ['Sp0246+1625', 'Sp0020 3305']
    eso.set(NIGHT, 'Europa', 100, '2026-10-09T13:55:00.000Z')
    poll(cfg, store, eso, at(14))
    eso.set(NIGHT, 'Europa', 300, '2026-10-09T14:25:00.000Z')
    poll(cfg, store, eso, at(14, 30))
    assert p0_jobs(store) == []
    poll(cfg, store, eso, at(15))                            # same count and last-modified for 30 min
    jobs = p0_jobs(store)
    assert [j['kind'] for j in jobs] == ['download', 'pipeline']
    dl, pl = jobs
    assert pl['depends_on'] == dl['id']
    assert 'sem:eso' in dl['locks'] and 'night:Europa:20261008' in dl['locks'] and not dl['disk_heavy']
    assert set(pl['locks']) == {'tel:Europa', 'night:Europa:20261008', 'target:SP0246+1625', 'target:SP0020--3305'}
    assert pl['argv'][:2] == ['./main/ZLP_pipeline.sh', '--force-platesolve']
    assert pl['argv'][2:] == ['1', cfg['basedir'], NIGHT, '8', '2', 'Europa']
    assert dl['argv'][-6:] == ['jobs', 'download', '--telescope', 'Europa', '--night', NIGHT]
    w = store.watch_get('Europa', NIGHT)
    assert w['state'] == 'enqueued' and 'unchanged' in w['ready_reason']
    before = len(eso.calls)
    poll(cfg, store, eso, at(15, 30))
    assert len(p0_jobs(store)) == 2                         # queued once
    assert ('summary', NIGHT) not in eso.calls[before:]      # a queued night is not asked about again


def test_transfer_log_match_queues_at_once(cfg, store):
    cfg = only(cfg, 'Ganymede')
    transfer_log(cfg, 'Ganymede', [('20261007', 136), (NIGHT, 390)])
    eso = FakeEso()
    eso.set(NIGHT, 'Ganymede', 390, '2026-10-09T14:51:19.677Z')
    poll(cfg, store, eso, at(14, 55))
    assert [j['kind'] for j in p0_jobs(store)] == ['download', 'pipeline']
    assert 'transfer log 390' in store.watch_get('Ganymede', NIGHT)['ready_reason']


def test_below_the_transfer_log_count_it_keeps_waiting(cfg, store):
    cfg = only(cfg, 'Europa')
    transfer_log(cfg, 'Europa', [(NIGHT, 874)])
    eso = FakeEso()
    eso.set(NIGHT, 'Europa', 370, '2026-10-09T14:53:14.593Z')
    for t in (at(15), at(15, 30), at(16), at(17)):
        poll(cfg, store, eso, t)
    assert p0_jobs(store) == []
    assert store.watch_get('Europa', NIGHT)['state'] == 'arriving'
    poll(cfg, store, eso, at(18, 5))                         # unchanged for over 180 min: take what is there
    assert [j['kind'] for j in p0_jobs(store)] == ['download', 'pipeline']
    assert 'look-back fetches the rest' in store.watch_get('Europa', NIGHT)['ready_reason']


def test_no_frames_and_no_transfer_entry_is_one_quiet_note(cfg, store):
    cfg = only(cfg, 'Io')
    transfer_log(cfg, 'Io', [('20260829', 499)])
    eso = FakeEso()
    def notes():
        return [e for e in store.events(kinds=['note']) if e['night'] == NIGHT]
    for t in (at(12), at(15), at(18)):
        poll(cfg, store, eso, t)
    assert notes() == []                                      # still within the day frames could arrive
    for t in (at(21, 30), at(22), at(22, 30)):
        poll(cfg, store, eso, t)
    assert len(notes()) == 1 and 'nothing to download' in notes()[0]['message']
    assert len(store.events(kinds=['note'])) == 2             # the night before last got its own single note
    assert store.jobs() == []
    assert store.watch_get('Io', NIGHT)['state'] == 'no_data'


def test_frames_turning_up_after_the_note_are_still_queued(cfg, store):
    cfg = only(cfg, 'Io')
    eso = FakeEso()
    poll(cfg, store, eso, at(21, 30))
    eso.set(NIGHT, 'Io', 50, '2026-10-09T21:40:00Z')
    poll(cfg, store, eso, at(22, 30))
    poll(cfg, store, eso, at(23, 0))
    assert [j['kind'] for j in p0_jobs(store)] == ['download', 'pipeline']


def test_eso_failure_changes_nothing(cfg, store):
    cfg = only(cfg, 'Europa')
    eso = FakeEso()
    eso.fail = True
    out = poll(cfg, store, eso, at(14))
    assert store.jobs() == [] and store.watch_get('Europa', NIGHT) is None
    assert any('failed' in line for line in out)
    assert store.events(kinds=['eso-error'])


def test_artemis_manifest_queues_the_pipeline_without_a_download(cfg, store):
    cfg = only(cfg, 'Artemis')
    d = os.path.join(cfg['basedir'], 'Observations', 'Artemis', 'images', NIGHT)
    names = ['Sp0020+3305-S001-R001-C00{}-i.fts'.format(k) for k in range(1, 4)] + \
            ['Sp0449+5138-S001-R001-C001-i.fts']
    for n in names:
        write_fits(os.path.join(d, n), OBJECT=n.split('-S0')[0], IMAGETYP='Light Frame')
    write_fits(os.path.join(d, 'AutoFlat', 'AutoFlat-Dusk-I+z-Bin1-001.fts'), IMAGETYP='Flat Field')
    with open(os.path.join(d, 'Data_Download.txt'), 'w') as f:
        f.write('\n'.join(names + ['AutoFlat/AutoFlat-Dusk-I+z-Bin1-001.fts']) + '\n')
    old = (at(10, 45) - dt.datetime(1970, 1, 1)).total_seconds()
    for p in [os.path.join(d, 'Data_Download.txt'), d]:
        os.utime(p, (old, old))
    poll(cfg, store, FakeEso(), at(11, 30))
    jobs = p0_jobs(store, 'Artemis')
    assert [j['kind'] for j in jobs] == ['pipeline']
    assert jobs[0]['depends_on'] is None and '--force-platesolve' not in jobs[0]['argv']
    assert {'target:SP0020+3305', 'target:SP0449+5138'} <= set(jobs[0]['locks'])
    assert 'all present' in store.watch_get('Artemis', NIGHT)['ready_reason']


def test_artemis_without_a_directory_gets_one_note_at_the_deadline(cfg, store):
    cfg = only(cfg, 'Artemis')
    poll(cfg, store, FakeEso(), at(15))
    assert [e for e in store.events(kinds=['note']) if e['night'] == NIGHT] == []
    poll(cfg, store, FakeEso(), at(22, 10))
    poll(cfg, store, FakeEso(), at(22, 40))
    assert len([e for e in store.events(kinds=['note']) if e['night'] == NIGHT]) == 1


def test_a_night_the_cron_already_ran_is_not_queued_again(cfg, store):
    cfg = only(cfg, 'Ganymede')
    transfer_log(cfg, 'Ganymede', [('20261007', 136)])
    logs = os.path.join(cfg['basedir'], 'PipelineOutput', 'v2', 'Ganymede', 'logs')
    os.makedirs(logs)
    open(os.path.join(logs, '20261007_1_v3.log'), 'w').close()      # T12 moved the cron run's log to v2
    eso = FakeEso()
    eso.set('20261007', 'Ganymede', 136, '2026-10-08T14:00:00Z')
    poll(cfg, store, eso, at(14))
    assert store.jobs() == []
    assert store.watch_get('Ganymede', '20261007')['state'] == 'processed'
    before = len(eso.calls)
    poll(cfg, store, eso, at(14, 30))
    assert ('summary', '20261007') not in eso.calls[before:]
