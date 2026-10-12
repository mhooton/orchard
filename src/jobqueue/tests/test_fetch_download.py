"""The verified add-only fetch, night requests from `q add` and backlog files, and the ESO bulk semaphore."""
import csv
import gzip
import io
import os

import pytest

from jobqueue import cli, jobs, jobspec
from jobqueue.fetch import EsoFileFetcher, Manifest, add_only_fetch
from jobqueue.nights import (FitsError, ForeignMatcher, diff_eso_disk, fits_layout, foreign_frames, foreign_times,
                             is_cube_layout, verify_frame)
from jobqueue.scheduler import plan
from jobqueue.tests.fakes import card, write_fits



def image_bytes(nx=4, ny=3, value=7, **keywords):
    """A complete 16-bit FITS image with data."""
    cards = [card('SIMPLE', True), card('BITPIX', 16), card('NAXIS', 2), card('NAXIS1', nx), card('NAXIS2', ny)]
    cards += [card(k.replace('_', '-'), v) for k, v in keywords.items()]
    cards.append('END'.ljust(80))
    head = ''.join(cards)
    head += ' ' * (-len(head) % 2880)
    data = bytes([0, value]) * (nx * ny)
    data += b'\0' * (-len(data) % 2880)
    return head.encode('ascii') + data


def cube_bytes(n=2):
    cards = [card('SIMPLE', True), card('BITPIX', 16), card('NAXIS', 3), card('NAXIS1', 4), card('NAXIS2', 3),
             card('NAXIS3', n), card('EXTEND', True), 'END'.ljust(80)]
    head = ''.join(cards)
    head += ' ' * (-len(head) % 2880)
    data = b'\0\1' * (4 * 3 * n)
    data += b'\0' * (-len(data) % 2880)
    ext = [card('XTENSION', 'BINTABLE'), card('BITPIX', 8), card('NAXIS', 2), card('NAXIS1', 10),
           card('NAXIS2', n), card('PCOUNT', 0), card('GCOUNT', 1), card('TFIELDS', 0), card('EXTNAME', 'METADATA'),
           'END'.ljust(80)]
    ehead = ''.join(ext)
    ehead += ' ' * (-len(ehead) % 2880)
    edata = b'x' * (10 * n)
    edata += b'\0' * (-len(edata) % 2880)
    return head.encode('ascii') + data + ehead.encode('ascii') + edata


def put(path, blob):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'wb') as f:
        f.write(blob)
    return str(path)


# ---------------------------------------------------------------- FITS checks


def test_verify_accepts_complete_frames_and_cubes(tmp_path):
    frame = put(tmp_path / 'a.fits', image_bytes())
    cube = put(tmp_path / 'c.fits', cube_bytes())
    assert [h['naxis'] for h in verify_frame(frame)] == [[4, 3]]
    assert not is_cube_layout(verify_frame(frame))
    assert is_cube_layout(verify_frame(cube))


@pytest.mark.parametrize('blob, why', [
    (image_bytes()[:-100], 'multiple of 2880'),                       # truncated mid-block
    (image_bytes()[:2880], 'need'),                                   # the data block is missing
    (b'NOTFITS ' + image_bytes()[8:], 'SIMPLE'),
    (image_bytes() + b'garbage!' * 360, 'XTENSION'),                  # a block that is no HDU
])
def test_verify_rejects_broken_files(tmp_path, blob, why):
    with pytest.raises(FitsError, match=why):
        verify_frame(put(tmp_path / 'bad.fits', blob))


def test_verify_rejects_a_header_only_file(tmp_path):
    with pytest.raises(FitsError, match='primary image'):
        verify_frame(write_fits(str(tmp_path / 'h.fits'), OBJECT='x'))
    assert fits_layout(str(tmp_path / 'h.fits'))[0]['data_bytes'] == 0


# ---------------------------------------------------------------- ACP names


def test_acp_named_frames_count_as_present_by_date_obs(tmp_path):
    night = tmp_path / '20180905'
    put(night / 'Sp1609-3431-S001-R001-C001-I+z.fts', image_bytes(DATE_OBS='2018-09-05T23:11:08.710'))
    put(night / 'Calibration' / 'Bias-S001-R001-C001-B1.fts', image_bytes(DATE_OBS='2018-09-06T10:25:04.480'))
    put(night / 'SPECU1.2018-09-05T23:30:00.000.fits', image_bytes())
    assert sorted(foreign_frames(str(night))) == ['Calibration/Bias-S001-R001-C001-B1.fts',
                                                  'Sp1609-3431-S001-R001-C001-I+z.fts']
    rows = [{'dp_id': 'SPECU1.2018-09-05T23:11:08.710'}, {'dp_id': 'SPECU1.2018-09-06T10:25:04.700'},
            {'dp_id': 'SPECU1.2018-09-05T23:30:00.000'}, {'dp_id': 'SPECU1.2018-09-05T23:40:00.000'}]
    present, missing, _ = diff_eso_disk(rows, os.listdir(night), foreign_times(str(night)))
    assert [r['dp_id'] for r in missing] == ['SPECU1.2018-09-05T23:40:00.000']
    assert len(present) == 3                         # exact, within 0.5 s, and by name


def test_a_foreign_frame_matches_one_row_only():
    m = ForeignMatcher({'2018-09-05T23:11:08.710': 'a.fts'})
    assert m.find('2018-09-05T23:11:08.710') == 'a.fts'
    assert m.find('2018-09-05T23:11:08.900') is None


# ---------------------------------------------------------------- downloading


class Response:
    def __init__(self, status=200, body=b'', headers=None, length=None):
        self.status_code = status
        self.body = body
        self.headers = dict(headers or {})
        if length is not None:
            self.headers['Content-Length'] = str(length)

    def iter_content(self, chunk_size=1):
        for k in range(0, len(self.body), chunk_size):
            yield self.body[k:k + chunk_size]

    def close(self):
        pass


class Token:
    def __init__(self):
        self.forced = 0

    def __call__(self, force=False):
        self.forced += bool(force)
        return 'tok'


def fetcher(responses, **kw):
    """responses: dp_id -> list of Response (one per attempt)."""
    calls = []

    def get(url, headers, timeout):
        dp_id = url.rsplit('/', 1)[1]
        calls.append(dp_id)
        assert headers['Authorization'] == 'Bearer tok'
        return responses[dp_id].pop(0)
    f = EsoFileFetcher(kw.pop('token', Token()), get=get, workers=2, attempts=3, backoff=0, log=lambda m: None,
                       sleep=lambda s: None, **kw)
    return f, calls


def test_fetcher_verifies_decompresses_and_retries(tmp_path):
    good = image_bytes()
    gz = io.BytesIO()
    with gzip.GzipFile(fileobj=gz, mode='wb') as z:
        z.write(good)
    f, calls = fetcher({
        'SPECU2.A': [Response(body=good, length=len(good))],
        'SPECU2.GZ': [Response(body=gz.getvalue(), headers={'Content-Type': 'application/gzip'})],
        'SPECU2.CUT': [Response(body=good[:3000], length=len(good)), Response(body=good)],   # cut, then whole
        'SPECU2.GONE': [Response(status=404)],
        'SPECU2.BROKEN': [Response(body=good[:-10]) for _ in range(3)],
    })
    got = f(['SPECU2.A', 'SPECU2.GZ', 'SPECU2.CUT', 'SPECU2.GONE', 'SPECU2.BROKEN'], str(tmp_path))
    assert sorted(os.path.basename(p) for p in got) == ['SPECU2.A.fits', 'SPECU2.CUT.fits', 'SPECU2.GZ.fits']
    for p in got:
        assert open(p, 'rb').read() == good
    assert calls.count('SPECU2.CUT') == 2 and calls.count('SPECU2.GONE') == 1 and calls.count('SPECU2.BROKEN') == 3
    assert set(f.failed) == {'SPECU2.GONE', 'SPECU2.BROKEN'}
    assert f.retry_later == set()                 # neither will work later either
    assert sorted(os.listdir(tmp_path)) == ['SPECU2.A.fits', 'SPECU2.CUT.fits', 'SPECU2.GZ.fits']   # no leftovers


def test_fetcher_renews_the_token_and_stops_after_repeated_network_failures(tmp_path):
    token = Token()
    down = [Response(status=503) for _ in range(3)]
    f, calls = fetcher({'SPECU2.H': [Response(headers={'Content-Type': 'text/html'}), Response(body=image_bytes())],
                        'SPECU2.X': list(down), 'SPECU2.Y': list(down), 'SPECU2.Z': list(down)},
                       token=token, breaker=4)
    f.workers = 1
    got = f(['SPECU2.H', 'SPECU2.X', 'SPECU2.Y', 'SPECU2.Z'], str(tmp_path))
    assert [os.path.basename(p) for p in got] == ['SPECU2.H.fits'] and token.forced == 1
    assert f.retry_later == {'SPECU2.X', 'SPECU2.Y', 'SPECU2.Z'}
    assert 'not tried' in f.failed['SPECU2.Z'] or calls.count('SPECU2.Z') < 3


def test_a_token_failure_is_retried_later_not_a_crash(tmp_path):
    def token(force=False):
        raise RuntimeError('ESO authentication failed')
    f, calls = fetcher({}, token=token)
    assert f(['SPECU2.A'], str(tmp_path)) == [] and f.retry_later == {'SPECU2.A'} and calls == []
    assert 'no ESO token' in f.failed['SPECU2.A']


# ---------------------------------------------------------------- placing: repair, refetch, manifest


def stage_download(blobs):
    def download(dp_ids, dest):
        return [put(os.path.join(dest, d + '.fits'), blobs[d]) for d in dp_ids if d in blobs]
    return download


def test_repair_swaps_only_frames_that_fail_verification(tmp_path):
    night = tmp_path / 'images' / '20250412'
    good_old = image_bytes(value=1)
    put(night / 'SPECU2.T.fits', image_bytes(value=2)[:-500])       # truncated
    put(night / 'SPECU2.OK.fits', good_old)
    manifest = tmp_path / 'manifest.csv'
    res = add_only_fetch([{'dp_id': 'SPECU2.T', 'dp_type': 'OBJECT'}, {'dp_id': 'SPECU2.OK'},
                          {'dp_id': 'SPECU2.NEW', 'dp_type': 'BIAS'}],
                         str(night), str(tmp_path / 'st'),
                         stage_download({d: image_bytes(value=3) for d in ('SPECU2.T', 'SPECU2.OK', 'SPECU2.NEW')}),
                         verify=verify_frame, replace='invalid', quarantine=str(tmp_path / 'q'),
                         record=Manifest(str(manifest), telescope='Europa', night='20250412', job='9'),
                         log=lambda m: None)
    assert res['replaced'] == 1 and res['added'] == 1 and res['already_present'] == 1
    assert open(night / 'SPECU2.T.fits', 'rb').read() == image_bytes(value=3)
    assert open(night / 'SPECU2.OK.fits', 'rb').read() == good_old                   # valid: left alone
    assert open(tmp_path / 'q' / 'SPECU2.T.fits', 'rb').read() == image_bytes(value=2)[:-500]   # kept
    rows = list(csv.DictReader(open(manifest)))
    assert sorted((r['action'], r['file']) for r in rows) == [('added', 'SPECU2.NEW.fits'),
                                                              ('replaced', 'SPECU2.T.fits')]
    assert all(r['md5'] and r['telescope'] == 'Europa' and r['job'] == '9' for r in rows)
    assert not (tmp_path / 'st').exists() and sorted(os.listdir(night)) == [
        'SPECU2.NEW.fits', 'SPECU2.OK.fits', 'SPECU2.T.fits']


def test_refetch_swaps_frames_that_differ_and_keeps_identical_ones(tmp_path):
    night = tmp_path / 'n'
    put(night / 'SPECU4.SAME.fits', image_bytes(value=5))
    put(night / 'SPECU4.DIFF.fits', image_bytes(value=6))
    put(night / 'SPECU4.FZ.fits.fz', b'compressed')
    res = add_only_fetch([{'dp_id': d} for d in ('SPECU4.SAME', 'SPECU4.DIFF', 'SPECU4.FZ')], str(night),
                         str(tmp_path / 'st'), stage_download({d: image_bytes(value=5) for d in
                                                               ('SPECU4.SAME', 'SPECU4.DIFF', 'SPECU4.FZ')}),
                         verify=verify_frame, replace='all', quarantine=str(tmp_path / 'q'), log=lambda m: None)
    assert res['replaced'] == 1 and res['already_present'] == 2 and res['added'] == 0
    assert os.listdir(tmp_path / 'q') == ['SPECU4.DIFF.fits']
    assert not (night / 'SPECU4.FZ.fits').exists()                     # never a second copy next to the .fz


def test_add_only_never_duplicates_an_acp_named_frame(tmp_path):
    night = tmp_path / 'n'
    put(night / 'Sp1609-3431-S001-R001-C001-I+z.fts', image_bytes(DATE_OBS='2018-09-05T23:11:08.710'))
    blob = image_bytes(DATE_OBS='2018-09-05T23:11:08.710')
    res = add_only_fetch([{'dp_id': 'SPECU1.2018-09-05T23:11:08.710'}], str(night), str(tmp_path / 'st'),
                         stage_download({'SPECU1.2018-09-05T23:11:08.710': blob}), verify=verify_frame,
                         foreign=foreign_times(str(night)), log=lambda m: None)
    assert res['added'] == 0 and res['foreign_duplicates'] == 1
    assert os.listdir(night) == ['Sp1609-3431-S001-R001-C001-I+z.fts']


def test_a_frame_that_fails_verification_is_not_placed(tmp_path):
    def broken(path):
        raise FitsError('nope')
    res = add_only_fetch([{'dp_id': 'SPECU2.A'}], str(tmp_path / 'n'), str(tmp_path / 'st'),
                         stage_download({'SPECU2.A': image_bytes()}), verify=broken, log=lambda m: None)
    assert res['invalid'] == ['SPECU2.A.fits'] and res['settled'] == [] and not (tmp_path / 'n').exists()


def test_replace_needs_a_quarantine_directory(tmp_path):
    with pytest.raises(ValueError):
        add_only_fetch([], str(tmp_path), str(tmp_path / 'st'), stage_download({}), replace='invalid')


# ---------------------------------------------------------------- the fetch job without a rows list


class Archive:
    def __init__(self, rows):
        self.rows = rows
        self.asked = []

    def frames(self, prog, first, last):
        self.asked.append((prog, first, last))
        return [dict(r) for r in self.rows if first <= r.get('night', first) <= last]


def test_download_only_fetch_lists_the_night_itself_and_queues_no_rerun(cfg, store, monkeypatch):
    night = '20250412'
    rows = [{'dp_id': 'SPECU2.2025-04-12T23:00:0{}.000'.format(k), 'dp_type': 'OBJECT', 'object': 'SP1234-5678',
             'origfile': '', 'access_estsize': '5000', 'night': night} for k in range(3)]
    d = os.path.join(cfg['basedir'], 'Observations', 'Europa', 'images', night)
    put(os.path.join(d, rows[0]['dp_id'] + '.fits'), image_bytes())
    archive = Archive(rows)
    asked = []
    monkeypatch.setattr(jobs, '_archive', lambda cfg: archive)
    monkeypatch.setattr(jobs, 'fetch_rows', lambda cfg, arch, tel, n, todo, **kw: asked.append(todo) or {
        'requested': len(todo), 'downloaded': len(todo), 'added': len(todo), 'added_science': len(todo),
        'settled': [], 'failed': {}, 'retry_later': []})
    assert jobs.fetch(cfg, store, 'Europa', night, None, 'none') == 0
    assert archive.asked == [('60.A-9009(B)', night, night)]
    assert [r['dp_id'] for r in asked[0]] == [r['dp_id'] for r in rows[1:]]
    assert store.jobs() == []                                         # download only: nothing queued


def test_download_only_fetch_asks_to_be_retried_after_network_trouble(cfg, store, monkeypatch):
    monkeypatch.setattr(jobs, '_archive', lambda cfg: Archive([{'dp_id': 'SPECU2.A', 'dp_type': 'OBJECT'}]))
    monkeypatch.setattr(jobs, 'fetch_rows', lambda *a, **kw: {
        'requested': 1, 'downloaded': 0, 'added': 0, 'added_science': 0, 'settled': [],
        'failed': {'SPECU2.A': 'ConnectionError'}, 'retry_later': ['SPECU2.A']})
    assert jobs.fetch(cfg, store, 'Europa', '20250412', None, 'none') == jobs.EX_TEMPFAIL


# ---------------------------------------------------------------- q add: nights and backlog files


def test_parse_backlog_reads_the_manual_format_and_rejects_bad_lines(tmp_path, cfg):
    f = tmp_path / 'backlog.csv'
    f.write_text('\n'.join([
        'Callisto,20240815,1,0,1',
        'Ganymede,20250101,1,1,1,"Sp0000-0000 Sp1111-1111"',
        'Artemis,20260530,0,0,1',
        '# a comment',
        'Iooooo,20240101,0,0,1',
        'Artemis,20250521,0,1,Sp1428+3310',
        'Europa,20240101,0,1,1',
        'Artemis,20260531,1,0,1',
        'Europa,2024-01-01,1,0,1',
        'Io,20240101,0,0,0',
    ]))
    reqs, errors = cli.parse_backlog(str(f), cfg)
    assert [(r['telescope'], r['night'], r['download'], r['process'], r['replace']) for r in reqs] == [
        ('Callisto', '20240815', True, True, 'never'), ('Ganymede', '20250101', True, True, 'all'),
        ('Artemis', '20260530', False, True, 'never')]
    assert reqs[1]['run_targets'] == ['Sp0000-0000', 'Sp1111-1111']
    assert len(errors) == 6
    assert any('unknown telescope' in e for e in errors) and any('0 or 1' in e for e in errors)
    assert any('DELETE needs DOWNLOAD' in e for e in errors) and any('do not come from ESO' in e for e in errors)
    assert any('YYYYMMDD' in e for e in errors) and any('nothing to do' in e for e in errors)


class Eso:
    def __init__(self, rows):
        self.rows = rows

    def frames(self, prog, first, last):
        return [dict(r) for r in self.rows if first <= r['night'] <= last]


def eso_rows(night, n, obj='SP0246+1625'):
    return [{'dp_id': 'SPECU2.{}-{}-{}T23:00:{:02d}.000'.format(night[:4], night[4:6], night[6:], k),
             'dp_type': 'OBJECT', 'object': obj, 'origfile': '', 'access_estsize': '5000', 'night': night}
            for k in range(n)]


def test_a_download_request_queues_a_fetch_and_the_pipeline_after_it(cfg, store):
    nights = cli.EsoNights(cfg, Eso(eso_rows('20250412', 4)))
    req = {'telescope': 'Europa', 'night': '20250412', 'download': True, 'process': True, 'replace': 'never'}
    fetch_id = cli.queue_night(cfg, store, nights, req, 'P2', log=lambda m: None)
    fetch = store.get_job(fetch_id)
    assert fetch['kind'] == 'fetch' and fetch['class'] == 'P2'
    assert fetch['argv'][-2:] == ['--rerun-class', 'none'] and '--rows' not in fetch['argv']
    assert set(fetch['locks']) == {'night:Europa:20250412', 'sem:eso', 'sem:eso_bulk'}
    (pipe,) = store.jobs(kind='pipeline')
    assert pipe['depends_on'] == fetch_id and pipe['dep_requires_success'] == 1
    assert 'target:SP0246+1625' in pipe['locks']


def test_download_only_and_nothing_to_fetch(cfg, store):
    rows = eso_rows('20250412', 2)
    for r in rows:
        put(os.path.join(cfg['basedir'], 'Observations', 'Europa', 'images', '20250412', r['dp_id'] + '.fits'),
            image_bytes())
    nights = cli.EsoNights(cfg, Eso(rows + eso_rows('20250413', 2)))
    out = []
    cli.queue_night(cfg, store, nights, {'telescope': 'Europa', 'night': '20250412', 'download': True,
                                         'process': True, 'replace': 'never'}, 'P2', log=out.append)
    assert 'nothing to fetch' in out[0] and [j['kind'] for j in store.jobs()] == ['pipeline']
    cli.queue_night(cfg, store, nights, {'telescope': 'Europa', 'night': '20250413', 'download': True,
                                         'process': False, 'replace': 'never'}, 'P2', log=out.append)
    assert sorted(j['kind'] for j in store.jobs()) == ['fetch', 'pipeline']


def test_dry_run_queues_nothing(cfg, store):
    nights = cli.EsoNights(cfg, Eso(eso_rows('20250412', 3)))
    out = []
    cli.queue_night(cfg, store, nights, {'telescope': 'Europa', 'night': '20250412', 'download': True,
                                         'process': True, 'replace': 'never'}, 'P2', dry_run=True, log=out.append)
    assert store.jobs() == [] and 'would fetch' in out[0] and 'would run the pipeline' in out[0]


def test_from_file_queues_nothing_if_any_line_is_bad(cfg, store, tmp_path, monkeypatch, capsys):
    f = tmp_path / 'b.csv'
    f.write_text('Europa,20250412,1,0,1\nIooooo,20240101,0,0,1\n')
    monkeypatch.setattr(cli, '_eso', lambda cfg: pytest.fail('must not query ESO'))
    a = cli.build_parser().parse_args(['add', '--class', 'P2', '--from-file', str(f)])
    assert cli.cmd_add(cfg, store, a) == 2 and store.jobs() == []
    assert 'nothing queued' in capsys.readouterr().out


# ---------------------------------------------------------------- scheduling


class Snap:
    load5, mem_available_gb, disk_util, cpu_count = 1.0, 400.0, {'sda': 5.0, 'sdb': 5.0}, 128


def test_bulk_fetches_never_take_both_eso_slots(cfg, store, clock):
    a = jobspec.add(store, jobspec.fetch_job(cfg, 'P2', 'Europa', '20250101', rerun=False))[0]
    b = jobspec.add(store, jobspec.fetch_job(cfg, 'P2', 'Io', '20250102', rerun=False))[0]
    queued = store.jobs(states=['queued'])
    first = plan(queued, [], Snap(), cfg, clock())
    started = [d.job['id'] for d in first if d.start]
    assert started == [a]                                        # one bulk fetch at a time
    store.mark_running(a)
    p0 = jobspec.add(store, jobspec.download_job(cfg, 'Callisto', '20261008', 100))[0]
    second = plan(store.jobs(states=['queued']), store.jobs(states=['running']), Snap(), cfg, clock())
    assert {d.job['id']: d.start for d in second} == {p0: True, b: False}   # the nightly download still runs
