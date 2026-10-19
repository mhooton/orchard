"""
Unit tests for the per-night shutter detection (calibration/shutter_detect.py)
and how its record reaches the de-smear (smear.shutter_state, correct_frame,
T4) and the notifier's digest.

Covers in particular:
  - the flat ramp test naming the shutter state and readout rate of each flat,
    and refusing faint or nearly saturated flats
  - evening and morning flats that agree settling the whole night; only the
    night's OWN flats counting, never those borrowed from another night
  - one-sided flats needing the science to confirm them; science-only nights
  - a night whose flats disagree: the science bins locating the switch, and the
    bisection never placing it outside the bracket it returns
  - READOUTM narrowing the trail hypotheses and vetoing a ramp it contradicts
  - disagreements with CAMSHUTT, an override window and silent headers recorded,
    carried into SHUTD* cards and warned about
  - the measurement deciding the de-smear only for unambiguous segments and
    only with detect_auto_apply; the override table still first
  - the record written at T2, read back, refreshed; a failure never raising
  - the digest listing disagreements of nights reduced in the last 24 h

Synthetic frames are 640 x 1004 (600 x 1000 after trimming) with their own
instrument parameters, so a night is quick to build and measure.

Run:  pytest src/tests/test_shutter_detect.py -q
"""
import datetime as dt
import json
import os
import sys

import numpy as np
import pytest
from astropy.io import fits

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from calibration import pipeutils                   # noqa: E402
# imported before any test patches pipeutils: pipeflat binds get_instrument_parameters at import
from calibration import pipeflat                    # noqa: E402
from calibration import shutter_detect as sd        # noqa: E402
from calibration import smear                       # noqa: E402

TRIM = (600, 1000)
N_SHIFTS = 1000
SMALL = {
    'gain': 1.0, 'saturation_threshold': 62000, 'naxis2': 640,
    'overscan': {'top_row': -15, 'bottom_row': None, 'left_col': 4, 'right_col': None},
    'trim': {'top_row': 22, 'bottom_row': 622, 'left_col': 2, 'right_col': 1002},
    'smear': {'axis': 1, 'direction': 1, 'readout_time': {'NORMAL': 4.33, 'FAST': 1.5},
              'readout_rate_time': {'1MHz': 4.33, '3MHz': 1.5}, 'n_shifts': N_SHIFTS},
}
TEL = 'SPECULOOS-GANYMEDE'
NIGHT = '20260810'


@pytest.fixture
def small(monkeypatch):
    """Small frames with their own instrument parameters, small trail bins, no shipped overrides."""
    def params(hdul, trimmed=False):
        return SMALL
    monkeypatch.setattr(pipeutils, 'get_instrument_parameters', params)
    monkeypatch.setattr(sd, 'get_instrument_parameters', params)
    monkeypatch.setattr(sd, 'TRAIL_MIN_FRAMES', 6)
    monkeypatch.setattr(sd, 'TRAIL_SAMPLE', 6)
    monkeypatch.setattr(sd, 'TRAIL_BIN', 3)
    monkeypatch.setattr(sd, 'TRAIL_MIN_STARS', 10)
    monkeypatch.setattr(sd, 'SKY_SAMPLE', 6)
    monkeypatch.setattr(sd, 'SKY_MIN_FRAMES', 4)
    monkeypatch.setattr(smear, '_shipped_overrides', lambda: ())
    smear._warned.clear()
    yield
    smear._warned.clear()


_STARS = {}


def stars(sky=600.0):
    """Nine isolated Gaussian stars, one per 60 rows, inside the trail test's margins."""
    if sky not in _STARS:
        yy, xx = np.mgrid[:TRIM[0], :TRIM[1]]
        img = np.full(TRIM, sky)
        for i, y in enumerate(range(60, 541, 60)):
            img += 20000.0 * np.exp(-((yy - y) ** 2 + (xx - (330 + 37 * (i % 6))) ** 2) / (2 * 1.6 ** 2))
        _STARS[sky] = img
    return _STARS[sky]


def t_read_of(state):
    return {'shuttered': 0.0, '1MHz': 4.33, '3MHz': 1.5}[state]


def raw_frame(signal, rng, pedestal=320.0, noise=6.2):
    raw = np.full((640, 1004), pedestal)
    raw[22:622, 2:1002] += signal
    return np.round(raw + rng.normal(0, noise, raw.shape)).astype(np.float32)


def write(d, t, kind, exptime, signal, rng, camshutt=None, readoutm=None):
    os.makedirs(d, exist_ok=True)
    hdu = fits.PrimaryHDU(raw_frame(signal, rng))
    hdu.header['TELESCOP'] = TEL
    hdu.header['IMAGETYP'] = {'flat': 'Flat Field', 'science': 'Light Frame', 'bias': 'Bias Frame'}[kind]
    hdu.header['EXPTIME'] = exptime
    hdu.header['FILTER'] = 'I+z'
    hdu.header['DATE-OBS'] = t.strftime('%Y-%m-%dT%H:%M:%S.000')
    if camshutt is not None:
        hdu.header['CAMSHUTT'] = camshutt
    if readoutm is not None:
        hdu.header['READOUTM'] = readoutm
    hdu.header['JD'] = 2440587.5 + (t - dt.datetime(1970, 1, 1)).total_seconds() / 86400.0
    path = os.path.join(d, 'SPECU3.%s.fits' % hdu.header['DATE-OBS'])
    hdu.writeto(path, overwrite=True)
    return path


def flat_signal(state, exptime, level=20000.0):
    flat = np.full(TRIM, level)
    t = t_read_of(state)
    return smear.add_smear(flat, smear.smear_tau(t, N_SHIFTS, exptime)) if t else flat


def science_signal(state, exptime, sky=600.0):
    t = t_read_of(state)
    img = stars(sky)
    return smear.add_smear(img, smear.smear_tau(t, N_SHIFTS, exptime)) if t else img


def make_night(root, night=NIGHT, evening=None, morning=None, science=(), exptime=10.0, n_flats=3, seed=0,
               camshutt=None, readoutm=None, sky=600.0):
    """A raw night directory: evening flats at 23:00, science from 00:00 at 10-minute steps (one state per
    frame), morning flats at 10:00 the next day. States: 'shuttered', '1MHz', '3MHz'."""
    rng = np.random.default_rng(seed)
    d = os.path.join(str(root), 'images', night)
    day = dt.datetime.strptime(night, '%Y%m%d')
    for i in range(n_flats):
        if evening:
            e = 5.0 + i
            write(d, day + dt.timedelta(hours=23, seconds=30 * i), 'flat', e, flat_signal(evening, e), rng,
                  camshutt, readoutm)
        if morning:
            e = 5.0 + i
            write(d, day + dt.timedelta(days=1, hours=10, seconds=30 * i), 'flat', e, flat_signal(morning, e),
                  rng, camshutt, readoutm)
    for k, state in enumerate(science):
        write(d, day + dt.timedelta(days=1, minutes=10 * k), 'science', exptime,
              science_signal(state, exptime, sky), rng, camshutt, readoutm)
    return d


# ------------------------------------------------------------------ flats

@pytest.mark.parametrize('state, verdict, rate', [('shuttered', 'shuttered', None), ('1MHz', 'shutterless', '1MHz'),
                                                  ('3MHz', 'shutterless', '3MHz')])
def test_flat_ramp_names_state_and_rate(tmp_path, small, state, verdict, rate):
    rng = np.random.default_rng(1)
    path = write(str(tmp_path), dt.datetime(2026, 8, 10, 23), 'flat', 6.0, flat_signal(state, 6.0), rng)
    m = sd.measure_flat(path, TEL, NIGHT)
    assert m['verdict'] == verdict and m['rate'] == rate
    assert m['t_read'] == pytest.approx(t_read_of(state), abs=0.05)


@pytest.mark.parametrize('level, reason', [(800.0, 'too faint'), (55000.0, 'near saturation')])
def test_unusable_flats_are_not_judged(tmp_path, small, level, reason):
    rng = np.random.default_rng(2)
    path = write(str(tmp_path), dt.datetime(2026, 8, 10, 23), 'flat', 6.0, flat_signal('1MHz', 6.0, level), rng)
    m = sd.measure_flat(path, TEL, NIGHT)
    assert m['verdict'] == 'inconclusive' and m['reason'] == reason


def test_sides_split_at_the_science():
    sci = [{'t': '2026-08-11T00:00:00'}, {'t': '2026-08-11T08:00:00'}]
    flats = [{'t': '2026-08-10T23:00:00'}, {'t': '2026-08-11T00:30:00'}, {'t': '2026-08-11T07:50:00'},
             {'t': '2026-08-11T10:00:00'}]
    sd.assign_sides(flats, sci)
    assert [f['side'] for f in flats] == ['evening', 'evening', 'morning', 'morning']


def test_sides_without_science_split_at_the_largest_gap():
    flats = [{'t': '2026-08-10T23:00:00', 'night': NIGHT}, {'t': '2026-08-10T23:05:00', 'night': NIGHT},
             {'t': '2026-08-11T10:00:00', 'night': NIGHT}]
    sd.assign_sides(flats, [])
    assert [f['side'] for f in flats] == ['evening', 'evening', 'morning']
    one = [{'t': '2026-08-11T10:00:00', 'night': NIGHT}, {'t': '2026-08-11T10:05:00', 'night': NIGHT}]
    sd.assign_sides(one, [])
    assert {f['side'] for f in one} == {'morning'}


def _flat(verdict, t_read=None, exptime=6.0):
    t_read = {'shutterless': 4.3, 'shuttered': 0.06}.get(verdict) if t_read is None else t_read
    return {'verdict': verdict, 'rate': None, 't_read': t_read, 'exptime': exptime, 't': 'x'}


def test_a_twilight_needs_every_conclusive_flat_to_agree():
    assert sd.side_summary([_flat('shutterless'), _flat('shutterless'), _flat('inconclusive')])['state'] \
        == 'shutterless'
    assert sd.side_summary([_flat('shutterless'), _flat('shuttered')])['state'] == 'mixed'
    assert sd.side_summary([_flat('shuttered')])['state'] is None          # one flat is not enough
    s = sd.side_summary([_flat('shutterless', 4.3), _flat('shutterless', 1.5), _flat('shutterless', 4.3, 10.0)])
    assert s['state'] == 'shutterless' and s['rate'] is None and 'not one readout rate' in s['rate_note']


def test_twilight_rate_from_the_t_read_exposure_line():
    """Ganymede 2026-02-20, morning: the twilight adds ~0.018 s per s of exposure; the intercept is 1 MHz.
    Per flat, the 35 s and 60 s ones are outside 15% of 4.33 s."""
    morning = [(8, 4.59), (8, 4.59), (11, 4.66), (11, 4.65), (15, 4.75), (15, 4.75), (22, 4.87), (22, 4.89),
               (35, 5.10), (35, 5.12), (60, 5.56), (60, 5.54)]
    s = sd.side_summary([_flat('shutterless', t, e) for e, t in morning])
    assert s['rate'] == '1MHz' and s['t_read_intercept'] == pytest.approx(4.45, abs=0.05)
    assert s['slope'] == pytest.approx(0.0185, abs=0.002)
    # shuttered flats: the flat field's own structure, +0.009 s per s, from ~0
    shut = [(6, 0.055), (10, 0.09), (15, 0.13), (25, 0.21), (40, 0.355)]
    assert sd.ramp_intercept([_flat('shuttered', t, e) for e, t in shut])[0] == pytest.approx(0.0, abs=0.02)
    # an early morning twilight on a shuttered night (Ganymede 20260308): a sky gradient, not a readout
    dawn = [(27, 2.50), (27, 3.17), (37, 5.33), (37, 6.43), (37, 7.71)]
    s = sd.side_summary([_flat('shutterless', t, e) for e, t in dawn])
    assert s['rate'] is None and 'sky gradient' in s['rate_note']
    # a cloud spoils one flat: the robust line barely moves
    cloudy = morning[:4] + [(15, 4.20)] + morning[5:]
    assert sd.side_summary([_flat('shutterless', t, e) for e, t in cloudy])['rate'] == '1MHz'


# ------------------------------------------------------------------ the night

def test_agreeing_twilights_settle_the_night(tmp_path, small):
    d = make_night(tmp_path, evening='1MHz', morning='1MHz', science=['1MHz'] * 3)
    rec = sd.detect_night(d, NIGHT, TEL)
    r = rec['result']
    assert (r['state'], r['source'], r['confidence'], r['rate'], r['t_read']) == \
        ('shutterless', 'flats', 'unambiguous', '1MHz', 4.33)
    assert rec['science'] == []                              # the flats settle it: no science measured
    [seg] = rec['segments']
    assert seg['start'] == '2026-08-10T23:00:00' and seg['end'] == '2026-08-11T10:01:00'
    assert rec['disagreements'][0]['kind'] == 'no-header'     # shutterless and no header says so


def test_shuttered_night(tmp_path, small):
    d = make_night(tmp_path, evening='shuttered', morning='shuttered', science=['shuttered'] * 3,
                   camshutt=True)
    rec = sd.detect_night(d, NIGHT, TEL)
    assert rec['result']['state'] == 'shuttered' and rec['result']['confidence'] == 'unambiguous'
    assert rec['result']['t_read'] is None and rec['disagreements'] == []


def test_borrowed_flats_are_not_the_nights_own(tmp_path, small):
    # 20260809 has shuttered flats; 20260810 none of its own and shutterless science (Ganymede 20250408's case)
    make_night(tmp_path, '20260809', evening='shuttered', morning='shuttered')
    d = make_night(tmp_path, NIGHT, science=['1MHz'] * 6)
    rec = sd.detect_night(d, NIGHT, TEL)
    assert rec['counts'].get('flat') is None
    assert (rec['result']['state'], rec['result']['source'], rec['result']['rate']) == \
        ('shutterless', 'science', '1MHz')


def test_one_twilight_needs_the_science(tmp_path, small):
    rec = sd.detect_night(make_night(tmp_path, evening='1MHz', science=['1MHz'] * 6), NIGHT, TEL)
    assert (rec['result']['state'], rec['result']['source'], rec['result']['confidence']) == \
        ('shutterless', 'flats+science', 'unambiguous')
    # long exposures: no trail test, and no shuttered reference night for the sky ramp
    rec = sd.detect_night(make_night(tmp_path / 'b', evening='1MHz', science=['1MHz'] * 6, exptime=60.0),
                          NIGHT, TEL)
    assert (rec['result']['state'], rec['result']['confidence']) == ('shutterless', 'probable')


def test_one_twilight_contradicted_by_the_science(tmp_path, small):
    rec = sd.detect_night(make_night(tmp_path, evening='shuttered', science=['1MHz'] * 6), NIGHT, TEL)
    assert rec['result']['state'] is None and rec['result']['confidence'] == 'inconclusive'
    # the flats themselves are still known
    [seg] = rec['segments']
    assert seg['state'] == 'shuttered' and seg['end'] < '2026-08-11'


def test_science_alone_is_never_unambiguous_by_default(tmp_path, small, monkeypatch):
    d = make_night(tmp_path, science=['1MHz'] * 6)
    assert sd.detect_night(d, NIGHT, TEL)['result']['confidence'] == 'probable'
    monkeypatch.setattr(sd, 'SCIENCE_ALONE_SIGMA', 4.0)          # an opt-in margin (this night: ~4.6)
    assert sd.detect_night(d, NIGHT, TEL)['result']['confidence'] == 'unambiguous'


def test_trail_cannot_tell_shuttered_from_3mhz_without_readoutm(tmp_path, small, monkeypatch):
    monkeypatch.setattr(sd, 'TRAIL_SYS_FLOOR', 0.8)
    d = make_night(tmp_path, science=['shuttered'] * 6)
    assert sd.detect_night(d, NIGHT, TEL)['result']['confidence'] == 'inconclusive'
    d = make_night(tmp_path / 'b', science=['shuttered'] * 6, readoutm=2)      # Astra: 1 MHz
    rec = sd.detect_night(d, NIGHT, TEL)
    assert rec['result']['state'] == 'shuttered' and rec['science'][0]['hypotheses'] == ['1MHz']


def test_split_night_is_located_by_the_science(tmp_path, small):
    science = ['shuttered'] * 14 + ['1MHz'] * 16
    d = make_night(tmp_path, evening='shuttered', morning='1MHz', science=science)
    rec = sd.detect_night(d, NIGHT, TEL)
    r = rec['result']
    assert (r['state'], r['source'], r['confidence']) == ('split', 'flats+science', 'unambiguous')
    first_shutterless = (dt.datetime(2026, 8, 11) + dt.timedelta(minutes=10 * 14)).isoformat()
    last_shuttered = (dt.datetime(2026, 8, 11) + dt.timedelta(minutes=10 * 13)).isoformat()
    assert r['switch'][0] <= last_shuttered and r['switch'][1] >= first_shutterless
    assert sd.segment_for(rec, '2026-08-10T23:00:00.000')['state'] == 'shuttered'
    assert sd.segment_for(rec, r['switch'][0])['state'] == 'shuttered'
    assert sd.segment_for(rec, r['switch'][1])['state'] == 'shutterless'
    assert sd.segment_for(rec, '2026-08-11T10:00:00')['state'] == 'shutterless'


def test_one_twilight_and_a_switch_in_the_science(tmp_path, small):
    # evening flats at 1 MHz, no morning flats, the shutter back in after 9 frames: the whole-night stack
    # alone would read an average, the end bins show the switch
    d = make_night(tmp_path, evening='1MHz', science=['1MHz'] * 9 + ['shuttered'] * 9)
    rec = sd.detect_night(d, NIGHT, TEL)
    r = rec['result']
    assert (r['state'], r['confidence']) == ('split', 'probable')
    early = sd.segment_for(rec, '2026-08-10T23:00:00')
    late = sd.segment_for(rec, '2026-08-11T02:50:00')
    assert (early['state'], early['confidence'], early['t_read']) == ('shutterless', 'unambiguous', 4.33)
    # the late bins rule out 1 MHz but cannot tell shuttered from 3 MHz within the PSF allowance
    assert (late['state'], late['confidence']) == (None, 'inconclusive') and 'rules out 1MHz' in late['note']


def test_switch_found_when_the_shuttered_end_cannot_exclude_3mhz(tmp_path, small, monkeypatch):
    # the PSF allowance hides 3 MHz from 0 s, so the shuttered bins are inconclusive; they still rule out
    # the 1 MHz the late bins read as, and that is enough to find and bracket the switch
    monkeypatch.setattr(sd, 'TRAIL_SYS_FLOOR', 0.6)
    d = make_night(tmp_path, evening='shuttered', science=['shuttered'] * 9 + ['1MHz'] * 9)
    rec = sd.detect_night(d, NIGHT, TEL)
    r = rec['result']
    assert r['state'] == 'split' and r['switch'][1] >= '2026-08-11T01:30:00'
    early = sd.segment_for(rec, '2026-08-10T23:00:00')
    late = sd.segment_for(rec, '2026-08-11T02:50:00')
    assert (early['state'], early['confidence']) == ('shuttered', 'unambiguous')
    assert (late['state'], late['rate'], late['confidence']) == ('shutterless', '1MHz', 'probable')
    # an allowance so wide that not even the 1 MHz end is conclusive: never "unambiguous shuttered"
    monkeypatch.setattr(sd, 'TRAIL_SYS_FLOOR', 0.8)
    rec = sd.detect_night(d, NIGHT, TEL)
    assert rec['result']['confidence'] != 'unambiguous'


def test_a_switch_on_a_night_without_flats(tmp_path, small):
    d = make_night(tmp_path, science=['shuttered'] * 9 + ['1MHz'] * 9)
    rec = sd.detect_night(d, NIGHT, TEL)
    r = rec['result']
    assert (r['state'], r['source'], r['confidence']) == ('split', 'science', 'probable')
    assert r['switch'][0] <= '2026-08-11T01:20:00' and r['switch'][1] >= '2026-08-11T01:30:00'
    # never a conclusive whole-night verdict from the mixture
    assert all(s['confidence'] != 'unambiguous' for s in rec['segments'])


@pytest.mark.parametrize('switch', range(0, 61, 3))
def test_bisection_always_brackets_the_switch(monkeypatch, switch):
    """Whatever a bin straddling the switch reads as, the switch stays inside the bracket."""
    frames = [{'t': str(i), 'i': i, 'exptime': 10.0} for i in range(60)]
    monkeypatch.setattr(sd, 'TRAIL_BIN', 6)
    for straddle in ('first', 'last', 'neither'):
        def stack(fs, telescope, night, log):
            after = [f['i'] >= switch for f in fs]
            if any(after) and not all(after) and straddle != 'last':
                if straddle == 'neither':
                    return {'verdict': 'inconclusive', 'rate': None, 'z': {'shuttered': 4.0, '1MHz': 2.0}}
                return {'verdict': 'inconclusive', 'rate': None, 'z': {'shuttered': 1.0, '1MHz': 9.0}}
            if any(after):
                return {'verdict': 'shutterless', 'rate': '1MHz', 'z': {'shuttered': 9.0, '1MHz': 0.0}}
            # shuttered: the trail test cannot rule out 3 MHz, but it rules out 1 MHz
            return {'verdict': 'inconclusive', 'rate': None, 'z': {'shuttered': 0.5, '1MHz': 9.0}}
        monkeypatch.setattr(sd, 'trail_stack', stack)
        _, found = sd.locate_switch(frames, TEL, NIGHT)
        if switch < 6 or switch > 54:
            continue                                         # an end bin straddles: nothing to bisect
        assert found is not None, (straddle, switch)
        lo, hi = found['lo'], found['hi']
        assert found['early'] is None and found['late'] == ('shutterless', '1MHz')
        assert lo < switch <= hi, (straddle, switch, found)
        assert hi - lo <= 2 * 6 + 6


def test_readoutm_contradicting_the_ramp_is_flagged(tmp_path, small):
    d = make_night(tmp_path, evening='1MHz', morning='1MHz', readoutm=1, camshutt=False)  # Astra says 3 MHz
    rec = sd.detect_night(d, NIGHT, TEL)
    assert rec['result']['confidence'] == 'probable'
    [dis] = [x for x in rec['disagreements'] if x['kind'] == 'READOUTM']
    assert dis['measured'] == '1MHz' and dis['header'] == ['3MHz']


def test_camshutt_disagreement_is_recorded_and_carded(tmp_path, small, capsys):
    d = make_night(tmp_path, evening='1MHz', morning='1MHz', science=['1MHz'] * 3, camshutt=True)
    rec = sd.detect_night(d, NIGHT, TEL)
    [dis] = rec['disagreements']
    assert dis['kind'] == 'CAMSHUTT' and dis['header'] == 'T' and dis['n_frames'] == dis['of'] == 9
    sd.warn_disagreements(rec)
    assert 'CAMSHUTT = T' in capsys.readouterr().out
    h = fits.Header({'DATE-OBS': '2026-08-11T03:00:00.000', 'CAMSHUTT': True, 'TELESCOP': TEL})
    cards = {k: v for k, v, _ in sd.header_cards(rec, h)}
    assert cards == {'SHUTDET': 'SHUTTERLESS', 'SHUTDSRC': 'flats', 'SHUTDCNF': 'unambiguous',
                     'SHUTDTRD': 4.33, 'SHUTDIS': True}


def test_override_disagreement_is_recorded(tmp_path, small):
    d = make_night(tmp_path, evening='shuttered', morning='shuttered')
    rows = [{'telescope': TEL, 'start': dt.datetime(2026, 8, 10, 12), 'end': dt.datetime(2026, 8, 11, 12),
             'camshutt': False, 'readout_time': None, 'note': ''}]
    rec = sd.detect_night(d, NIGHT, TEL, overrides=rows)
    [dis] = rec['disagreements']
    assert dis['kind'] == 'override' and dis['override'] == 'F' and dis['measured'] == 'shuttered'


def test_nothing_usable_is_inconclusive(tmp_path, small):
    d = make_night(tmp_path, science=['1MHz'] * 3)              # too few frames for the trail test
    rec = sd.detect_night(d, NIGHT, TEL)
    assert rec['result']['confidence'] == 'inconclusive' and rec['segments'] == []
    h = fits.Header({'DATE-OBS': '2026-08-11T03:00:00.000'})
    assert {k: v for k, v, _ in sd.header_cards(rec, h)}['SHUTDET'] == 'UNKNOWN'


def test_sky_ramp_against_a_shuttered_reference_night(tmp_path, small):
    make_night(tmp_path, '20260809', evening='shuttered', morning='shuttered', n_flats=3)
    d = make_night(tmp_path, NIGHT, science=['1MHz'] * 6, exptime=60.0, sky=3000.0)
    rec = sd.detect_night(d, NIGHT, TEL)
    [sky] = rec['science']
    assert sky['method'] == 'sky_ramp' and sky['reference_night'] == '20260809'
    assert sky['t_read'] == pytest.approx(4.33, abs=0.3)
    assert (rec['result']['state'], rec['result']['source'], rec['result']['confidence']) == \
        ('shutterless', 'sky_ramp', 'probable')                 # experimental: never unambiguous


def test_sky_ramp_only_when_nothing_else_answers(tmp_path, small):
    make_night(tmp_path, '20260809', evening='shuttered', morning='shuttered', n_flats=3)
    d = make_night(tmp_path, NIGHT, evening='1MHz', science=['1MHz'] * 6, exptime=60.0, sky=3000.0)
    assert sd.detect_night(d, NIGHT, TEL)['science'] == []             # the evening flats answered
    rec = sd.detect_night(d, NIGHT, TEL, all_levels=True)                # validation runs everything
    assert [s['method'] for s in rec['science']] == ['sky_ramp']


def test_judge():
    rates = {'1MHz': 4.33, '3MHz': 1.5}
    assert sd.judge(4.4, 0.1, 0.4, rates)[:2] == ('shutterless', '1MHz')
    assert sd.judge(1.45, 0.05, 0.2, rates)[:2] == ('shutterless', '3MHz')
    assert sd.judge(0.1, 0.05, 0.2, rates)[:2] == ('shuttered', None)
    assert sd.judge(0.1, 0.05, 0.8, rates)[0] == 'inconclusive'          # 3 MHz not excluded
    assert sd.judge(0.1, 0.05, 0.8, {'1MHz': 4.33})[0] == 'shuttered'
    assert sd.judge(9.0, 0.1, 0.4, rates)[0] == 'inconclusive'           # consistent with nothing
    assert sd.judge(float('nan'), 0.1, 0.4, rates)[0] == 'inconclusive'


def test_candidate_rates_follow_readoutm():
    rates = {'1MHz': 4.33, '3MHz': 1.5}
    assert sd.candidate_rates([{'readoutm': 2}, {'readoutm': 2}], rates) == ({'1MHz': 4.33}, '1MHz')
    assert sd.candidate_rates([{'readoutm': 2}, {'readoutm': None}], rates) == (rates, None)
    assert sd.candidate_rates([{'readoutm': 0}], rates) == ({}, '5MHz')    # no measured time at 5 MHz


# ------------------------------------------------------------------ records

def test_record_written_read_back_and_refreshed(tmp_path, small):
    d = make_night(tmp_path, evening='shuttered', morning='shuttered')
    reddir = str(tmp_path / 'output' / NIGHT / 'reduction')
    os.makedirs(reddir)
    rec = sd.load_or_detect(reddir, d, NIGHT, TEL)
    assert os.path.exists(os.path.join(reddir, 'shutter.json'))
    with open(os.path.join(reddir, 'shutter.json')) as f:
        assert json.load(f)['result'] == rec['result']
    rec2 = sd.load_or_detect(reddir, d, NIGHT, TEL)
    assert rec2['measured_at'] == rec['measured_at']
    # a borrowed night is measured on its flats only, into its own file
    sd.load_or_detect(reddir, d, '20260809', TEL, own=False)
    assert os.path.exists(os.path.join(reddir, 'shutter_20260809.json'))


def test_failed_measurement_never_raises(tmp_path, small):
    reddir = str(tmp_path)
    rec = sd.load_or_detect(reddir, None, NIGHT, TEL)
    assert rec['segments'] == [] and 'measurement failed' in rec['notes'][0]


def test_own_night_dir_from_the_lists(tmp_path):
    reddir = tmp_path / 'output' / NIGHT / 'reduction'
    reddir.mkdir(parents=True)
    (reddir / '1_flat_I+z.list').write_text('/x/images/20260809/a.fits\n/x/images/%s/AutoFlat/b.fits\n' % NIGHT)
    assert sd.own_night_dir(str(reddir), '1', NIGHT) == '/x/images/%s' % NIGHT
    assert sd.night_dir_of('/x/images/20260809/a.fits') == ('/x/images/20260809', '20260809')
    assert sd.night_dir_of('/x/elsewhere/a.fits') is None


def test_night_records_pick_the_frames_own_night(tmp_path, small):
    make_night(tmp_path, '20260809', evening='shuttered', morning='shuttered')
    own = make_night(tmp_path, NIGHT, evening='1MHz', morning='1MHz')
    reddir = str(tmp_path / 'output' / NIGHT / 'reduction')
    os.makedirs(reddir)
    records = sd.NightRecords(reddir, '1', refresh_borrowed=True)
    borrowed = os.path.join(str(tmp_path), 'images', '20260809', os.listdir(str(tmp_path / 'images' / '20260809'))[0])
    assert records.for_frame(borrowed, TEL)['result']['state'] == 'shuttered'
    assert records.for_frame(os.path.join(own, os.listdir(own)[0]), TEL)['result']['state'] == 'shutterless'


# ------------------------------------------------------------------ the de-smear decision

def _record(state='shutterless', confidence='unambiguous', t_read=4.33):
    return {'night': NIGHT, 'segments': [sd._segment('2026-08-10T23:00:00', '2026-08-11T10:00:00', state,
                                                     '1MHz' if t_read else None, t_read, 'flats', confidence)]}


def _hdr(camshutt=None, date_obs='2026-08-11T03:00:00.000'):
    h = fits.Header({'TELESCOP': TEL, 'DATE-OBS': date_obs, 'EXPTIME': 20.0})
    if camshutt is not None:
        h['CAMSHUTT'] = camshutt
    return h


def test_measurement_decides_only_with_auto_apply():
    rec = _record()
    assert smear.shutter_state(_hdr(True), [], rec, auto_apply=False) == (True, 'CAMSHUTT', None)
    assert smear.shutter_state(_hdr(True), [], rec, auto_apply=True) == (False, 'measured (flats)', 4.33)
    assert smear.shutter_state(_hdr(), [], _record('shuttered', t_read=None), True) == \
        (True, 'measured (flats)', None)


def test_measurement_needs_an_unambiguous_segment_and_loses_to_overrides():
    assert smear.shutter_state(_hdr(True), [], _record(confidence='probable'), True) == (True, 'CAMSHUTT', None)
    assert smear.shutter_state(_hdr(), [], _record(t_read=None), True) == (True, 'default', None)
    rows = [{'telescope': TEL, 'start': dt.datetime(2026, 8, 10, 12), 'end': dt.datetime(2026, 8, 11, 12),
             'camshutt': True, 'readout_time': None, 'note': ''}]
    assert smear.shutter_state(_hdr(), rows, _record(), True) == (True, 'override', None)
    # a frame outside every segment (another night's) falls through to CAMSHUTT
    assert smear.shutter_state(_hdr(False, '2026-08-12T03:00:00'), [], _record(), True) == \
        (False, 'CAMSHUTT', None)


def test_correct_frame_cards_and_measured_readout_time():
    andor = {'saturation_threshold': 62000, 'smear': dict(SMALL['smear'], n_shifts=2048)}
    img = stars()
    h = _hdr(True)
    smear._warned.clear()
    out = smear.correct_frame(img, h, andor, overrides=[], detection=_record())
    assert out is img and h['SMEARCOR'] is False and h['SMEARSRC'] == 'CAMSHUTT'
    assert h['SHUTDET'] == 'SHUTTERLESS' and h['SHUTDIS'] is True
    andor['smear']['detect_auto_apply'] = True
    h = _hdr(True)
    smear.correct_frame(img, h, andor, overrides=[], detection=_record())
    assert h['SMEARCOR'] is True and h['SMEARTRD'] == 4.33 and h['SMEARSRC'] == 'measured (flats)'
    assert '(measured)' in h.comments['SMEARTRD']


def test_shipped_config_does_not_auto_apply():
    cfg = pipeutils.load_instrument_config()
    for tel, instruments in cfg['telescopes'].items():
        for name, params in instruments.items():
            assert not smear.auto_apply(params), (tel, name)


def test_t4_desmears_flats_on_the_measurement_when_enabled(tmp_path, small, monkeypatch):
    params = dict(SMALL, smear=dict(SMALL['smear'], detect_auto_apply=True))
    monkeypatch.setattr(pipeutils, 'get_instrument_parameters', lambda hdul, trimmed=False: params)
    monkeypatch.setattr(sd, 'get_instrument_parameters', lambda hdul, trimmed=False: params)
    monkeypatch.setattr(pipeflat, 'get_instrument_parameters', lambda hdul, trimmed=False: params)
    d = make_night(tmp_path, evening='1MHz', morning='1MHz', n_flats=3)      # no CAMSHUTT, no override
    reddir = os.path.join(str(tmp_path), 'output', NIGHT, 'reduction') + '/'
    os.makedirs(reddir)
    fits.PrimaryHDU(np.zeros(TRIM, dtype=np.float32)).writeto(reddir + 'run_MasterBias.fits')
    fits.PrimaryHDU(np.zeros(TRIM, dtype=np.float32)).writeto(reddir + 'run_MasterDark.fits')
    flist = reddir + 'run_flat_I+z.list'
    with open(flist, 'w') as f:
        f.write('\n'.join(sorted(os.path.join(d, n) for n in os.listdir(d))) + '\n')
    pipeflat.reducer([flist], 'run_MasterBias.fits', 'run_MasterDark.fits', 'run_MasterFlat.fits', '1', '1',
                     reddir, reddir)
    with fits.open(reddir + 'run_MasterFlat_I+z.fits') as h:
        master, hdr = np.array(h[0].data), h[0].header
    assert hdr['NSMEARED'] == 6
    left, right = np.median(master[:, 100:300]), np.median(master[:, -300:-100])
    assert abs(right / left - 1) < 5e-3
    assert os.path.exists(reddir + 'shutter.json')


# ------------------------------------------------------------------ the digest

def test_digest_lists_recent_disagreements(tmp_path, monkeypatch):
    from jobqueue import digest as D
    monkeypatch.setattr(D, 'render_work', lambda q, now: [])

    class Ctx:
        basedir = str(tmp_path)
        now = dt.datetime(2026, 10, 18, 7, 45)

    def record(tel, night, disagreements, age_h):
        d = tmp_path / 'PipelineOutput' / 'v3' / tel / 'output' / night / 'reduction'
        d.mkdir(parents=True)
        p = d / 'shutter.json'
        p.write_text(json.dumps({'night': night, 'disagreements': disagreements}))
        t = (Ctx.now - dt.timedelta(hours=age_h) - dt.datetime(1970, 1, 1)).total_seconds()
        os.utime(str(p), (t, t))

    record('Ganymede', '20261016', [{'kind': 'CAMSHUTT', 'measured': 'shutterless', 'confidence': 'unambiguous',
                                     'header': 'T', 'n_frames': 812, 'of': 900}], 3)
    record('Io', '20261016', [{'kind': 'no-header', 'measured': 'shutterless', 'measured_t_read': 4.33,
                               'confidence': 'probable'}], 5)
    record('Europa', '20261010', [{'kind': 'CAMSHUTT', 'measured': 'shuttered', 'header': 'F', 'n_frames': 1,
                                   'of': 1}], 72)                      # too old for this digest
    lines = D.shutter_notes(Ctx)
    assert len(lines) == 2
    assert 'Ganymede, night of 16 Oct 2026: frames measure shutterless [unambiguous] but CAMSHUTT = T on 812' \
        in lines[0]
    assert 'Io' in lines[1] and '4.33 s' in lines[1] and 'not de-smeared' in lines[1]
    cfg = {'notify': {'subject_prefix': '[orchard]'}}
    q = {'failures': {}, 'lookback_actions': [], 'eso_errors': [], 'hours': 24}
    _, body = D.render_digest(cfg, '20261016', [], q, [], [], [], Ctx.now, shutter=lines)
    assert 'Shutter state measured against the headers' in body and lines[0] in body
