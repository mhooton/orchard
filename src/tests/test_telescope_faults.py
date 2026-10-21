"""
Tests for utils/telescope_faults.py: the frames a night lost to the telescope (a stopped mount, a pointing that
walked away), from raw headers, processed frames, the WCS cache and T7's log lines. The synthetic nights copy the
measured cases: Ganymede 2025-07-16 TOI-4336 c (Astra: TRACKING 0 from 00:30 UT, RA at 15 deg/h, airmass fixed,
CRVAL 0) and Artemis 2024-02-11 Sp1356+4342 (MaxIm DL: the pointing jumped 4.3' at 01:34 UT, no solution after
01:49). Headers are written with the standard library, so these run without astropy; the WCS arithmetic is also
checked against astropy where it is installed.
"""
import datetime as dt
import math
import os
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from utils import telescope_faults as tf  # noqa: E402

SCALE = 0.35 / 3600.0           # deg per pixel, Andor
NIGHT = '20250716'
T0 = dt.datetime(2025, 7, 16, 23, 0, 0)


# ---------------------------------------------------------------------------
# synthetic FITS
# ---------------------------------------------------------------------------

def _card(key, value):
    if isinstance(value, bool):
        text = '{:<8}= {:>20}'.format(key, 'T' if value else 'F')
    elif isinstance(value, (int, float)):
        text = '{:<8}= {:>20}'.format(key, repr(value))
    else:
        text = "{:<8}= '{:<8}'".format(key, str(value).replace("'", "''"))
    return text.ljust(80)[:80]


def _header_bytes(cards):
    text = ''.join(_card(k, v) for k, v in cards) + 'END'.ljust(80)
    return text.ljust(-(-len(text) // 2880) * 2880).encode('ascii')


def write_fits(path, cards, extension=None):
    with open(path, 'wb') as f:
        f.write(_header_bytes(cards))
        if extension is not None:
            f.write(_header_bytes(extension))


def wcs_cards(ra0, dec0, x, y, rot=0.0, parity=-1.0, scale=SCALE):
    """A TAN WCS putting (ra0, dec0) at 1-based pixel (x, y), the detector turned by `rot` degrees."""
    c, s = math.cos(math.radians(rot)), math.sin(math.radians(rot))
    cd = [[parity * scale * c, scale * s], [-parity * scale * s, scale * c]]
    return [('CTYPE1', 'RA---TAN'), ('CTYPE2', 'DEC--TAN'), ('CRVAL1', ra0), ('CRVAL2', dec0), ('CRPIX1', x),
            ('CRPIX2', y), ('CD1_1', cd[0][0]), ('CD1_2', cd[0][1]), ('CD2_1', cd[1][0]), ('CD2_2', cd[1][1])]


def astra_frame(d, i, t, obj='TOI-4336 c', ra=206.10513, dec=-40.3364, tracking=1.0, slewing=0.0, parked=0,
                airmass=None, wcs=True, kind='Light Frame'):
    name = 'SPECU3.{}.fits'.format(t.strftime('%Y-%m-%dT%H:%M:%S.000'))
    cards = [('SIMPLE', True), ('BITPIX', 16), ('NAXIS', 2), ('NAXIS1', 2048), ('NAXIS2', 2088),
             ('OBJECT', obj), ('IMAGETYP', kind), ('EXPTIME', 10), ('DATE-OBS', t.strftime('%Y-%m-%dT%H:%M:%S.000')),
             ('RA', ra), ('DEC', dec), ('AIRMASS', airmass if airmass is not None else 1.04 + 0.0004 * i),
             ('ALTITUDE', 63.3), ('AZIMUTH', 227.3), ('TRACKING', tracking), ('SLEWING', slewing),
             ('TELPARK', parked), ('ASTRA', '0.3.0')]
    if wcs:
        cards += wcs_cards(ra, dec, 752.0, 1038.0, rot=0.3)
    else:                                       # Astra's empty solution on a trailed frame
        cards += [('CTYPE1', 'nan'), ('CTYPE2', 'nan'), ('CRVAL1', 0.0), ('CRVAL2', 0.0), ('CRPIX1', 0.0),
                  ('CRPIX2', 2089.0)]
    write_fits(os.path.join(d, name), cards)
    return name


def ganymede_night(d, n_ok=40, n_stop=30, cadence=17):
    """TOI-4336 c: tracking, then the mount stops (TRACKING 0, RA running at the sidereal rate, airmass fixed), then
    another target with tracking resumed; a flat and a dark either side."""
    names = []
    astra_frame(d, 0, T0 - dt.timedelta(minutes=30), obj='', kind='FLAT', tracking=0.54, slewing=0.46, wcs=False)
    t = T0
    for i in range(n_ok):
        # the first frame starts while the mount still settles: TRACKING 0.74 + SLEWING 0.26 is not a stop
        tr, sl = (0.74, 0.26) if i == 0 else (1.0, 0.0)
        names.append(astra_frame(d, i, t, tracking=tr, slewing=sl))
        t += dt.timedelta(seconds=cadence)
    t_stop, ra = t, 206.10513
    for i in range(n_stop):
        ra_i = ra + tf.SIDEREAL_DEG_PER_H * (t - t_stop).total_seconds() / 3600.0
        names.append(astra_frame(d, n_ok + i, t, ra=ra_i, tracking=0.0, airmass=1.1187153, wcs=False))
        t += dt.timedelta(seconds=cadence)
    for i in range(20):
        astra_frame(d, i, t + dt.timedelta(minutes=10, seconds=69 * i), obj='Sp2135-3153', ra=323.81104,
                    dec=-31.89607)
    astra_frame(d, 0, t + dt.timedelta(hours=3), obj='', kind='Dark Frame', tracking=0.0, parked=1, wcs=False)
    return names, t_stop


def maxim_frame(d, i, t, obj='Sp1356+4342', ra='13 56 41.28', dec='+43 42 59.2', ext='fts'):
    name = '{}-S001-R001-C{:03d}-I+z.{}'.format(obj, i + 1, ext)
    write_fits(os.path.join(d, name), [
        ('SIMPLE', True), ('BITPIX', 16), ('NAXIS', 2), ('NAXIS1', 2048), ('NAXIS2', 2088), ('OBJECT', obj),
        ('IMAGETYP', 'Light Frame'), ('EXPTIME', 23.0), ('DATE-OBS', t.strftime('%Y-%m-%dT%H:%M:%S.380')),
        ('RA', ra), ('DEC', dec), ('OBJCTRA', ra), ('OBJCTDEC', dec), ('AIRMASS', 1.9),
        ('ALTITUDE', 30.8), ('OBJCTALT', 30.8), ('SWCREATE', 'MaxIm DL Version 6.14')])
    return name


RA_ART, DEC_ART = 15 * (13 + 56 / 60. + 41.28 / 3600.), 43 + 42 / 60. + 59.2 / 3600.


def artemis_night(d, procdir, n=60, jump_at=30, solved_until=40, wrong=(), drift_per_frame=0.0, jump_px=(-690, 265)):
    """Sp1356+4342: the pointing sits at pixel (1052, 1017) of the processed frames, jumps by jump_px at frame
    jump_at (or walks drift_per_frame px a frame from there), and nothing solves after solved_until. Frames in
    `wrong` carry a wrong old solution (turned 116 deg)."""
    os.makedirs(procdir, exist_ok=True)
    t = dt.datetime(2024, 2, 12, 0, 43, 59)
    names = []
    for i in range(n):
        names.append(maxim_frame(d, i, t))
        key = names[-1][:-4]
        x, y = 1052.0, 1017.0
        if i >= jump_at:
            x, y = x + jump_px[0] + drift_per_frame * (i - jump_at), y + jump_px[1]
        cards = [('SIMPLE', True), ('BITPIX', -32), ('NAXIS', 2), ('NAXIS1', 2046), ('NAXIS2', 2044),
                 ('OBJECT', 'Sp1356+4342')]
        if i in wrong:
            cards += wcs_cards(RA_ART + 0.05, DEC_ART - 0.02, 700.0, 1500.0, rot=116.4, parity=1.0)
        elif i < solved_until:
            cards += wcs_cards(RA_ART, DEC_ART, x, y, rot=0.6, parity=1.0)
        write_fits(os.path.join(procdir, 'proc' + key + '.fits'), cards)
        t += dt.timedelta(seconds=34)
    return names


def t7_line(name, ok=True, error='Plate scale sanity check failed: solved=10.6 arcsec/px', crpix=(751.9, 1015.9),
            crval=(323.791675, -31.892262), n=1, total=100):
    if ok:
        return '[{}/{}] proc{} SUCCESS - Sources: 332→16 Gaia: 100→32 Matches: 16 Method: pairs ' \
               'CRPIX=({:.2f},{:.2f}) CRVAL=({:.6f},{:.6f})'.format(n, total, name, crpix[0], crpix[1], crval[0],
                                                                     crval[1])
    return '[{}/{}] proc{} FAILED - Sources: 294→16 Gaia: 353→32 Matches: 0 Error: {}'.format(n, total, name, error)


# ---------------------------------------------------------------------------
# headers
# ---------------------------------------------------------------------------

def test_header_values(tmp_path):
    p = str(tmp_path / 'a.fits')
    write_fits(p, [('SIMPLE', True), ('NAXIS', 2), ('OBJECT', "TOI-4336 c"), ('NOTE', "it's"), ('EXPTIME', 10),
                   ('AIRMASS', 1.1187153272511539), ('TRACKING', 0.0), ('FLAG', False)])
    h = tf.read_header(p)
    assert h['OBJECT'] == 'TOI-4336 c' and h['NOTE'] == "it's"
    assert h['EXPTIME'] == 10 and h['AIRMASS'] == pytest.approx(1.1187153272511539) and h['TRACKING'] == 0.0
    assert h['SIMPLE'] is True and h['FLAG'] is False


def test_compressed_frame_reads_the_image_header(tmp_path):
    p = str(tmp_path / 'procX.fits.fz')
    write_fits(p, [('SIMPLE', True), ('BITPIX', 8), ('NAXIS', 0), ('EXTEND', True)],
               extension=[('XTENSION', 'BINTABLE'), ('NAXIS', 2), ('NAXIS1', 8), ('NAXIS2', 2044), ('ZIMAGE', True),
                          ('ZNAXIS1', 2046), ('ZNAXIS2', 2044)] + wcs_cards(10.0, -20.0, 1000.0, 1000.0))
    h = tf.read_header(p)
    assert h['NAXIS1'] == 2046 and h['NAXIS2'] == 2044 and h['CRVAL1'] == 10.0
    assert tf.header_wcs(h)['shape'] == (2044, 2046)


def test_angles_and_times():
    assert tf.angle('13 56 41.28', hours=True) == pytest.approx(209.172)
    assert tf.angle('+43 42 59.2') == pytest.approx(43.716444, abs=1e-6)
    assert tf.angle('-00 30 00') == pytest.approx(-0.5)
    assert tf.angle(206.10513) == 206.10513
    assert tf.frame_time({'DATE-OBS': '2026-09-14T22:55:43.072000'}) == dt.datetime(2026, 9, 14, 22, 55, 43, 72000)
    assert tf.frame_time({'DATE-OBS': '2024-02-12', 'TIME-OBS': '00:43:59'}) == dt.datetime(2024, 2, 12, 0, 43, 59)
    assert tf.frame_key('procSPECU3.2025-07-17T02:00:50.000.fits') == 'SPECU3.2025-07-17T02:00:50.000'
    assert tf.frame_key('Sp1356+4342-S001-R001-C001-I+z.fts') == 'Sp1356+4342-S001-R001-C001-I+z'
    assert tf.frame_key('procSp1356+4342-S001-R001-C001-I+z.fits.fz') == 'Sp1356+4342-S001-R001-C001-I+z'


def test_astra_empty_and_placeholder_wcs_are_none():
    assert tf.header_wcs({'CTYPE1': 'nan', 'CTYPE2': 'nan', 'CRVAL1': 0.0, 'CRVAL2': 0.0}) is None
    assert tf.header_wcs({'CTYPE1': 'PIXEL', 'CTYPE2': 'PIXEL', 'CRVAL1': 1.0, 'CRVAL2': 1.0}) is None
    pc = {'CTYPE1': 'RA---TAN', 'CTYPE2': 'DEC--TAN', 'CRVAL1': 10.0, 'CRVAL2': 5.0, 'CRPIX1': 1.0, 'CRPIX2': 1.0,
          'CDELT1': -SCALE, 'CDELT2': SCALE}
    assert tf.header_wcs(pc)['cd'] == [[-SCALE, 0.0], [0.0, SCALE]]


# ---------------------------------------------------------------------------
# the TAN projection
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('ra0,dec0,rot,parity', [(206.1, -40.3, 0.3, -1.0), (209.17, 43.72, 116.0, 1.0),
                                                  (359.99, 2.0, 270.0, -1.0), (45.0, -89.5, 10.0, 1.0)])
def test_projection_round_trip_and_astropy(ra0, dec0, rot, parity):
    w = tf.header_wcs(dict(wcs_cards(ra0, dec0, 1023.5, 1022.0, rot=rot, parity=parity)))
    for dx, dy in ((0, 0), (900, -700), (-1000, 1000)):
        ra, dec = tf.pixel_to_world(w, 1023.5 + dx, 1022.0 + dy)
        x, y = tf.world_to_pixel(w, ra, dec)
        assert (x, y) == pytest.approx((1023.5 + dx, 1022.0 + dy), abs=1e-6)
    assert tf.scale_arcsec(w) == pytest.approx(0.35)
    wcs_mod = pytest.importorskip('astropy.wcs')
    fits = pytest.importorskip('astropy.io.fits')
    h = fits.Header()
    for k, v in wcs_cards(ra0, dec0, 1023.5, 1022.0, rot=rot, parity=parity):
        h[k] = v
    aw = wcs_mod.WCS(h)
    for px, py in ((10.0, 2000.0), (1500.0, 300.0)):
        ra, dec = aw.all_pix2world([[px, py]], 1)[0]
        assert tf.pixel_to_world(w, px, py) == pytest.approx((ra, dec), abs=1e-9)
        assert tf.world_to_pixel(w, ra, dec) == pytest.approx((px, py), abs=1e-6)


def test_flip_turns_the_rotation_by_half_a_turn():
    a = tf.header_wcs(dict(wcs_cards(10.0, 5.0, 1.0, 1.0, rot=0.6)))
    b = tf.header_wcs(dict(wcs_cards(10.0, 5.0, 1.0, 1.0, rot=180.6)))
    assert tf._turn(tf.rotation_deg(a), tf.rotation_deg(b)) == pytest.approx(180.0)
    assert tf._turn(tf.rotation_deg(a), tf.rotation_deg(b), 180.0) == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------
# tracking stopped
# ---------------------------------------------------------------------------

def test_tracking_stop_is_found_with_its_frames_and_times(tmp_path):
    d = str(tmp_path)
    names, t_stop = ganymede_night(d)
    res = tf.night(d)
    assert set(res) == {'TOI-4336 c', 'Sp2135-3153'}
    toi = res['TOI-4336 c']
    assert toi['frames'] == 70 and toi['tracking_known']
    (e,) = toi['episodes']
    assert e['kind'] == 'tracking' and (e['i0'], e['i1'], e['n']) == (41, 70, 30)
    assert e['first'] == names[40] and e['last'] == names[-1]
    assert e['start'] == t_stop.isoformat()
    assert e['detail']['ra_rate'] == pytest.approx(tf.SIDEREAL_DEG_PER_H, abs=0.05)
    assert e['detail']['airmass_span'] < 1e-6 and e['detail']['parked'] == 0
    assert res['Sp2135-3153']['episodes'] == []
    assert tf.flags(res) == ('Y', 'N')
    text = tf.describe(e, toi['frames'])
    assert text.startswith('TOI-4336 c: tracking stopped 23:11–23:19 UT, frames 41–70 of 70 (30 frames)')
    assert 'header RA ran 15.0 deg/h, airmass fixed at 1.119' in text
    assert tf.digest_text(res) == 'telescope: tracking stopped on TOI-4336 c 23:11–23:19 UT (30 frames)'


def test_settling_frames_and_calibrations_are_not_stops(tmp_path):
    d = str(tmp_path)
    ganymede_night(d, n_stop=0)
    res = tf.night(d)
    assert tf.episodes(res) == [] and tf.flags(res) == ('N', 'N') and tf.digest_text(res) == ''


def test_parked_light_frames_are_stops(tmp_path):
    d = str(tmp_path)
    t = T0
    for i in range(6):
        astra_frame(d, i, t, tracking=0.0 if i >= 3 else 1.0, parked=1 if i >= 3 else 0, wcs=i < 3)
        t += dt.timedelta(minutes=1)
    (e,) = tf.night(d)['TOI-4336 c']['episodes']
    assert e['kind'] == 'tracking' and e['n'] == 3 and e['detail']['parked'] == 3
    assert '3 with the telescope parked' in tf.describe(e) and 'altitude' not in tf.describe(e)
    e['detail']['altitude'] = 20.4             # Ganymede 2025-06-21 Ch_SN2024ggi_z stopped at airmass 2.9
    assert tf.describe(e).endswith('at altitude 20 deg')


def test_tracking_zero_while_the_airmass_moves_is_no_stop(tmp_path):
    """Callisto 2024-12-09 Sp0903+0540: TRACKING 0 for half an hour while the airmass rose 2.031 -> 2.019, as only a
    tracking mount's does. Astra's flag was wrong; a stopped mount's airmass stays fixed."""
    d = str(tmp_path)
    t = T0
    for i in range(30):
        astra_frame(d, i, t, tracking=0.0 if i >= 10 else 1.0, airmass=2.031 - 0.0004 * i)
        t += dt.timedelta(minutes=1)
    assert tf.episodes(tf.night(d)) == []


def test_maxim_frames_cannot_say_but_trailed_solves_can(tmp_path):
    """Without TRACKING the headers cannot tell; T7's 'Sources trailed' failures can, three in a row."""
    d = str(tmp_path)
    t = dt.datetime(2024, 2, 12, 0, 43, 59)
    names = [maxim_frame(d, i, t + dt.timedelta(seconds=34 * i)) for i in range(10)]
    assert tf.flags(tf.night(d)) == ('unknown', 'unknown')
    log = tmp_path / 'log.txt'
    trail = 'Sources trailed (axis ratio 4.2, all one way), the telescope was probably not tracking; twirl: x'
    lines = [t7_line(n[:-4] + '.fits', ok=False, error=trail if 4 <= i <= 7 else 'Plate solve failed')
             for i, n in enumerate(names)]
    lines[5] = lines[5] + ' ' + lines[6]            # two workers printed on one line
    del lines[6]
    log.write_text('TARGET: Sp1356+4342\n' + '\n'.join(lines) + '\n')
    res = tf.night(d, [str(log)])
    (e,) = res['Sp1356+4342']['episodes']
    assert e['kind'] == 'trailed' and (e['i0'], e['i1']) == (5, 8)
    assert tf.flags(res)[0] == 'Y'
    assert 'stars trailed (T7)' in tf.describe(e)


# ---------------------------------------------------------------------------
# the mount's own position log (technical_logs/<YYYY-MM-DD>/<TEL>_pos.csv)
# ---------------------------------------------------------------------------

POS_HEAD = '#Telescope and Dome positions\n\n#jd,date,time,ALT,AZ,DEC,RA,DEC(inst.),HA(inst.),{}Dome_AZ\n'


def pos_log(path, start, end, stops=(), gaps=(), slews=(), target_dist=True, inst=True):
    """The mount's position log as the site writes it, a line a second from start to end. Tracking, the hour-angle
    axis turns at the sidereal rate and the altitude and azimuth follow; within each (t0, t1) of `stops` every axis
    stands still and the log's RA runs at 15 deg/h (Callisto 2024-10-12: altitude 58.6375, azimuth 36.5484 for
    2h20m); within each of `slews` the axes move a degree a second; no lines within `gaps`. `target_dist` False is
    the log before mid-2019; `inst` False gives 0 for both instrumental axes, as some nights' logs do."""
    rate = tf.SIDEREAL_DEG_PER_H / 3600.0
    alt, az, ha, ra = 58.2, 36.1, 18.0, 0.07
    lines = [POS_HEAD.format('Target_Dist,' if target_dist else '')]
    t = start
    while t <= end:
        if not any(a <= t < b for a, b in gaps):
            jd = 2440587.5 + (t - tf.EPOCH).total_seconds() / 86400.0
            cols = ['{:.8f}'.format(jd), t.strftime('%Y-%m-%d'), t.strftime('%H:%M:%S.000'), '{:.9f}'.format(alt),
                    '{:.9f}'.format(az), '1.260000000', '{:.9f}'.format(ra), '{:.9f}'.format(-1.26 if inst else 0.0),
                    '{:.9f}'.format(ha if inst else 0.0)] + (['0.296000005'] if target_dist else []) + \
                ['{:.9f}'.format(az)]
            lines.append(','.join(cols) + '\n')
        if any(a <= t < b for a, b in stops):
            ra += rate / 15.0
        elif any(a <= t < b for a, b in slews):
            ha, alt = ha - 1.0, alt + 0.5
        else:
            ha, alt, az = ha - rate, alt + 0.3 * rate, az + 0.6 * rate
        t += dt.timedelta(seconds=1)
    with open(path, 'w') as f:
        f.writelines(lines)


def acp_log(path, night_, slews, telescope='Europa', tracking_off=None):
    """ACP's log as the site delivers it: '<TEL>,<night's date>,<UT time>,<message>', the evening's date kept after
    midnight. A 'Start slew to <target>' line for each (time, target) of `slews`, by the ACP PC's clock."""
    day = '{}-{}-{}'.format(night_[:4], night_[4:6], night_[6:])
    lines = ['#TELESCOPE,DATE,TIME,MESSAGE\n']
    for t, target in slews:
        lines.append('{},{},{},  Start slew to {}...\n'.format(telescope, day, t.strftime('%H:%M:%S'), target))
        lines.append('{},{},{},  (slew complete)\n'.format(telescope, day,
                                                          (t + dt.timedelta(seconds=12)).strftime('%H:%M:%S')))
    if tracking_off:
        lines.append('{},{},{},**RUN STOPPED AT SPECIFIED QUIT TIME\n'.format(telescope, day,
                                                                             tracking_off.strftime('%H:%M:%S')))
        lines.append('{},{},{},  (turning tracking off for safety)\n'.format(telescope, day,
                                                                             tracking_off.strftime('%H:%M:%S')))
    with open(path, 'w') as f:
        f.writelines(lines)


def maxim_sso_frame(d, t, obj='Sp0002+0115', airmass=1.3, exptime=36.0):
    """A pre-Astra SSO light frame: MaxIm DL's header (no TRACKING; RA, ALTITUDE and AIRMASS are the target's) under
    ESO's name."""
    name = 'SPECU4.{}.fits'.format(t.strftime('%Y-%m-%dT%H:%M:%S.000'))
    write_fits(os.path.join(d, name), [
        ('SIMPLE', True), ('BITPIX', 16), ('NAXIS', 2), ('NAXIS1', 2048), ('NAXIS2', 2088), ('OBJECT', obj),
        ('IMAGETYP', 'Light Frame'), ('EXPTIME', exptime), ('DATE-OBS', t.strftime('%Y-%m-%dT%H:%M:%S.000')),
        ('RA', '00 02 06.05'), ('DEC', '+01 15 36.0'), ('AIRMASS', airmass), ('ALTITUDE', 50.3),
        ('SWCREATE', 'MaxIm DL Version 6.29')])
    return name


def maxim_stop_night(obsdir, n=60, stop=(20, 45), night_='20230612', telescope='Callisto', lead=None):
    """Laid out as ZLP has it: n MaxIm DL frames of Sp0002+0115 from 02:00 UT, one every 43 s (36 s exposures) with
    the target's airmass falling, in <obsdir>/images/<night>; the night's position log in technical_logs, with the
    mount standing still from 5 s before frame stop[0] (0-based) until 3 s after frame stop[1] - 1 ends. stop None:
    tracking all night. lead (s): the position log's clock runs this far ahead of the frames', and ACP's log beside
    it has the night's two slews, which the position log also has."""
    raw = os.path.join(obsdir, 'images', night_)
    os.makedirs(raw)
    t0 = dt.datetime(int(night_[:4]), int(night_[4:6]), int(night_[6:]), 2, 0, 0) + dt.timedelta(days=1)
    names = [maxim_sso_frame(raw, t0 + dt.timedelta(seconds=43 * i), airmass=1.3 - 0.002 * i) for i in range(n)]
    shift = dt.timedelta(seconds=lead or 0)
    stops = []
    if stop:
        stops = [(t0 + dt.timedelta(seconds=43 * stop[0] - 5) + shift,
                  t0 + dt.timedelta(seconds=43 * (stop[1] - 1) + 39) + shift)]
    log = tf.mount_log_path(obsdir, telescope, night_)
    os.makedirs(os.path.dirname(log))
    start, slews = t0 - dt.timedelta(minutes=5), []
    if lead is not None:
        acp = [(t0 - dt.timedelta(minutes=40), 'Sp2331-2749'), (t0 - dt.timedelta(minutes=4), 'Sp0002+0115')]
        acp_log(tf.acp_log_path(obsdir, telescope, night_), night_, acp, telescope=telescope)
        start, slews = t0 - dt.timedelta(minutes=45), [(a + shift, a + shift + dt.timedelta(seconds=20)) for a, _ in acp]
    pos_log(log, start + shift, t0 + dt.timedelta(seconds=43 * n + 300) + shift, stops=stops, slews=slews)
    return raw, log, names


def test_mount_log_finds_the_stop_maxim_headers_cannot_show(tmp_path):
    raw, log, names = maxim_stop_night(str(tmp_path))
    assert tf.flags(tf.night(raw)) == ('unknown', 'unknown')
    res = tf.night(raw, mount_log=log)
    r = res['Sp0002+0115']
    assert r['frames'] == 60 and r['tracking_known']
    (e,) = r['episodes']
    assert e['kind'] == 'tracking' and (e['i0'], e['i1'], e['n']) == (21, 45, 25)
    assert e['first'] == names[20] and e['last'] == names[44]
    d = e['detail']
    # where the mount stood: the log's altitude after 1,155 s of tracking at 0.3 of the sidereal rate
    assert d['source'] == 'mount log' and d['altitude'] == pytest.approx(59.6, abs=0.05) and d['parked'] == 0
    assert 'airmass' not in d and 'ra_rate' not in d          # MaxIm DL's are the target's, not the mount's
    assert tf.flags(res) == ('Y', 'unknown')
    assert tf.describe(e, r['frames']) == ("Sp0002+0115: tracking stopped 02:14–02:31 UT, frames 21–45 of 60 "
                                           "(25 frames); the mount's position log has it standing still at "
                                           "altitude 60 deg")
    assert tf.digest_text(res) == 'telescope: tracking stopped on Sp0002+0115 02:14–02:31 UT (25 frames)'


def test_mount_log_of_a_tracking_mount_says_no(tmp_path):
    raw, log, _ = maxim_stop_night(str(tmp_path), stop=None)
    res = tf.night(raw, mount_log=log)
    assert tf.episodes(res) == [] and tf.flags(res) == ('N', 'unknown')
    # a night without its log (Artemis after June 2025, SSO before 2019) cannot say, as before
    assert tf.flags(tf.night(raw, mount_log=str(tmp_path / 'missing_pos.csv'))) == ('unknown', 'unknown')


def test_mount_log_exposure_shares(tmp_path):
    p = str(tmp_path / 'Io_pos.csv')
    t = dt.datetime(2023, 8, 11, 2, 0, 0)
    minute = dt.timedelta(minutes=1)
    pos_log(p, t, t + 30 * minute, stops=[(t + 10 * minute, t + 20 * minute)], gaps=[(t + 25 * minute, t + 27 * minute)])
    log = tf.MountLog(p)
    assert log.exposure(t + 2 * minute, 60)['idle'] == 0.0
    still = log.exposure(t + 12 * minute, 60)
    assert still['idle'] == 1.0 and (still['alt'], still['az']) == pytest.approx((58.2 + 0.3 * 600 * 15.041 / 3600,
                                                                                  36.1 + 0.6 * 600 * 15.041 / 3600))
    assert log.exposure(t + 9.5 * minute, 60)['idle'] == pytest.approx(0.5, abs=0.02)   # stopped half way through
    assert log.exposure(t + 24 * minute + dt.timedelta(seconds=40), 60) is None    # the log covers 19 s of it
    assert log.exposure(t + 24 * minute, 60)['idle'] == 0.0                         # ... and 59 s of this one
    assert log.exposure(t - 5 * minute, 60) is None and log.exposure(t + 40 * minute, 60) is None


@pytest.mark.parametrize('target_dist,inst', [(False, True), (True, False)])
def test_mount_log_formats(tmp_path, target_dist, inst):
    """Before mid-2019 the log has no Target_Dist; some nights give 0 for both instrumental axes, which leaves the
    altitude and azimuth to say whether the mount moved."""
    p = str(tmp_path / 'Io_pos.csv')
    t = dt.datetime(2019, 2, 21, 2, 0, 0)
    minute = dt.timedelta(minutes=1)
    pos_log(p, t, t + 30 * minute, stops=[(t + 10 * minute, t + 20 * minute)], target_dist=target_dist, inst=inst)
    log = tf.MountLog(p)
    assert log.exposure(t + 2 * minute, 60)['idle'] == 0.0 and log.exposure(t + 12 * minute, 60)['idle'] == 1.0


def test_astra_headers_are_not_second_guessed_by_the_mount_log(tmp_path, monkeypatch):
    """Where Astra's TRACKING says, the log is not even read."""
    d = str(tmp_path)
    ganymede_night(d)
    monkeypatch.setattr(tf, 'read_mount_log', lambda path: pytest.fail('read for frames that have TRACKING'))
    (e,) = tf.night(d, mount_log=str(tmp_path / 'Ganymede_pos.csv'))['TOI-4336 c']['episodes']
    assert e['kind'] == 'tracking' and 'source' not in e['detail']


def test_mount_log_path():
    assert tf.mount_log_path('/data/Observations/Io', 'Io', '20230810') == \
        '/data/Observations/Io/technical_logs/2023-08-10/Io_pos.csv'
    assert tf.acp_log_path('/data/Observations/Io', 'Io', '20230810') == \
        '/data/Observations/Io/technical_logs/2023-08-10/Io_ACP.txt'


def test_clock_offset_from_the_slews_both_logs_record(tmp_path):
    """Europa 2019-05-05: ACP slewed to Sp1331-6513 at 23:37:03, to an offset at 23:38:15 and back at 23:39:18 by
    its clock, and to Sp1741+0940 at 04:30:12; the position log has each slew begin 31 s later by its own. A slew
    ACP logged that the position log lacks has no vote, and ACP's lines after midnight keep the evening's date."""
    pos, acp = str(tmp_path / 'Europa_pos.csv'), str(tmp_path / 'Europa_ACP.txt')
    t = dt.datetime(2019, 5, 5, 23, 30, 0)
    slews = [(t + dt.timedelta(seconds=423), 'Sp1331-6513'), (t + dt.timedelta(seconds=495), 'offset'),
             (t + dt.timedelta(seconds=558), 'Sp1331-6513'), (t + dt.timedelta(hours=5, seconds=12), 'Sp1741+0940'),
             (t + dt.timedelta(hours=5, minutes=30), 'Sp1741+0940')]
    log_slews = [(a + dt.timedelta(seconds=31), a + dt.timedelta(seconds=40)) for a, _ in slews[:4]]
    pos_log(pos, t, t + dt.timedelta(hours=5, minutes=40), slews=log_slews)
    acp_log(acp, '20190505', slews)
    log = tf.MountLog(pos)
    assert len(log.slews) == 4
    assert tf.acp_slews(acp)[3] == ((dt.datetime(2019, 5, 6, 4, 30, 12) - tf.EPOCH).total_seconds(), 'Sp1741+0940')
    assert log.clock_offset(acp) == pytest.approx(31.0, abs=0.01)
    acp_log(acp, '20190505', slews[:1])                     # one slew is not enough to move the clock
    assert log.clock_offset(acp) is None
    assert log.clock_offset(str(tmp_path / 'none_ACP.txt')) is None


def test_clock_offset_ignores_acps_pointing_corrections(tmp_path):
    """Io 2019-05-15: ACP slewed to Sp1327-3549 at 22:59:45 and again at 23:00:35 to correct its pointing, then to
    Sp1921-2915 at 03:38:25 and 03:39:19; the position log saw the acquisitions only, 1-2 s after ACP's lines. Each
    correction would otherwise pair with the acquisition before it, at -49 s, and tie with the true offset."""
    pos, acp = str(tmp_path / 'Io_pos.csv'), str(tmp_path / 'Io_ACP.txt')
    t = dt.datetime(2019, 5, 15, 22, 50, 0)
    slews = [(t + dt.timedelta(seconds=585), 'Sp1327-3549'), (t + dt.timedelta(seconds=635), 'Sp1327-3549'),
             (t + dt.timedelta(hours=4, seconds=2905), 'Sp1921-2915'),
             (t + dt.timedelta(hours=4, seconds=2959), 'Sp1921-2915')]
    pos_log(pos, t, t + dt.timedelta(hours=5), slews=[(slews[i][0] + dt.timedelta(seconds=2),
                                                        slews[i][0] + dt.timedelta(seconds=16)) for i in (0, 2)])
    acp_log(acp, '20190515', slews, telescope='Io')
    log = tf.MountLog(pos)
    assert log.clock_offset(acp) == pytest.approx(2.0, abs=0.01)
    # with the corrections under other names nothing tells the two offsets apart: no offset rather than a wrong one
    acp_log(acp, '20190515', [(a, 'target{}'.format(i)) for i, (a, _) in enumerate(slews)], telescope='Io')
    assert log.clock_offset(acp) is None


def test_clock_offset_keeps_the_frame_before_acp_turns_tracking_off(tmp_path):
    """Europa 2019-04-21: the last 30-s frame began at 10:08:40 and ACP turned tracking off at 10:09:20, both by
    ACP's clock; the position log, 33 s behind, has the mount still from 10:08:47. Read without the offset the
    frame looks stopped for 77% of it; with it, the mount tracked throughout."""
    raw = tmp_path / 'images' / '20190421'
    raw.mkdir(parents=True)
    t0 = dt.datetime(2019, 4, 22, 9, 55, 0)
    names = [maxim_sso_frame(str(raw), t0 + dt.timedelta(seconds=40 * i), obj='Sp1741+0940', exptime=30.0)
             for i in range(20)]
    off = t0 + dt.timedelta(seconds=40 * 19 + 30 + 10)       # ACP's clock: 10 s after the last frame ends
    slews = [(t0 - dt.timedelta(minutes=50), 'Sp1656-2046'), (t0 - dt.timedelta(minutes=5), 'Sp1741+0940')]
    lag = dt.timedelta(seconds=33)
    pos, acp = str(tmp_path / 'Europa_pos.csv'), str(tmp_path / 'Europa_ACP.txt')
    pos_log(pos, t0 - dt.timedelta(hours=1), t0 + dt.timedelta(minutes=30),
            slews=[(a - lag, a - lag + dt.timedelta(seconds=15)) for a, _ in slews],
            stops=[(off - lag, t0 + dt.timedelta(hours=1))])
    acp_log(acp, '20190421', slews, tracking_off=off)
    (e,) = tf.night(str(raw), mount_log=pos)['Sp1741+0940']['episodes']
    assert e['n'] == 1 and e['last'] == names[-1] and e['detail']['clock_offset'] is None
    res = tf.night(str(raw), mount_log=pos, acp_log=acp)
    assert tf.episodes(res) == [] and tf.flags(res)[0] == 'N'


def test_t7_outcomes_last_attempt_wins(tmp_path):
    a, b = tmp_path / 'old.log', tmp_path / 'new.log'
    a.write_text(t7_line('SPECU3.X.fits', ok=False) + '\n')
    b.write_text('noise\n' + t7_line('SPECU3.X.fits', ok=True, crpix=(10.0, 20.0)) + '\n')
    out = tf.t7_outcomes([str(a), str(b)])
    assert out['SPECU3.X']['state'] == 'SUCCESS' and out['SPECU3.X']['crpix'] == (10.0, 20.0)
    assert out['SPECU3.X']['gaia'] == 100 and out['SPECU3.X']['sources'] == 332
    assert tf.t7_outcomes([str(b), str(a)])['SPECU3.X']['state'] == 'FAILED'
    assert tf.t7_outcomes([str(tmp_path / 'missing.log')]) == {}


# ---------------------------------------------------------------------------
# pointing drift
# ---------------------------------------------------------------------------

def test_jump_is_found_from_processed_frames_despite_wrong_old_solutions(tmp_path):
    d, out = str(tmp_path / 'raw'), str(tmp_path / 'v3' / '20240211')
    os.makedirs(d)
    artemis_night(d, os.path.join(out, 'Sp1356+4342', '1'), wrong=(45, 52, 53, 58))
    processed = tf.processed_frames([out])
    assert len(processed) == 60
    res = tf.night(d, processed=processed)
    r = res['Sp1356+4342']
    assert not r['tracking_known'] and r['solved'] == 40
    (e,) = r['episodes']
    assert e['kind'] == 'drift' and (e['i0'], e['i1'], e['n']) == (31, 60, 30)
    det = e['detail']
    assert det['max_arcmin'] == pytest.approx(4.3, abs=0.1) and det['jump_arcmin'] == pytest.approx(4.3, abs=0.1)
    assert det['unsolved_after'] and 'off_field_at' not in det
    assert tf.flags(res) == ('unknown', 'Y')
    text = tf.describe(e, r['frames'])
    assert text == ("Sp1356+4342: pointing jumped 4.3' at 01:00 UT; 01:00–01:17 UT, frames 31–60 of 60 "
                    "(30 frames); no frame solved after 01:06 UT")
    assert tf.digest_text(res) == "telescope: pointing jumped 4.3' on Sp1356+4342 01:00–01:17 UT (30 frames)"


def test_walk_off_the_field_is_a_drift_not_a_jump(tmp_path):
    d, out = str(tmp_path / 'raw'), str(tmp_path / 'v3' / '20240211')
    os.makedirs(d)
    artemis_night(d, os.path.join(out, 'Sp1356+4342', '1'), n=80, jump_at=20, solved_until=80,
                  drift_per_frame=-40.0, jump_px=(0, 0))
    (e,) = tf.night(d, processed=tf.processed_frames([out]))['Sp1356+4342']['episodes']
    det = e['detail']
    assert e['kind'] == 'drift' and 'jump_arcmin' not in det and not det['unsolved_after']
    assert e['i0'] == 1 + 20 + math.ceil(1.0 * 60 / 0.35 / 40.0) and e['i1'] == 80
    assert det['off_field_at'] is not None and det['max_arcmin'] == pytest.approx(59 * 40 * 0.35 / 60, abs=0.1)
    text = tf.describe(e)
    assert "pointing drifted, up to 13.8' away by 01:28 UT" in text, text
    assert 'target off the field from' in text


def test_guiding_jitter_and_lone_bad_solutions_are_no_drift(tmp_path):
    d, out = str(tmp_path / 'raw'), str(tmp_path / 'v3' / '20240211')
    os.makedirs(d)
    artemis_night(d, os.path.join(out, 'Sp1356+4342', '1'), jump_at=10 ** 6, solved_until=60, wrong=(20, 40))
    res = tf.night(d, processed=tf.processed_frames([out]))
    assert tf.episodes(res) == [] and tf.flags(res) == ('unknown', 'N')


def track_night(d, positions, step_s=34, rot=0.6):
    """MaxIm frames of Sp1356+4342 with a raw WCS each, the target at the given pixel positions (None: unsolved)."""
    t = dt.datetime(2017, 12, 13, 0, 25, 45)
    for i, xy in enumerate(positions):
        name = maxim_frame(d, i, t + dt.timedelta(seconds=step_s * i))
        if xy is not None:
            cards = list(tf.read_header(os.path.join(d, name)).items()) + wcs_cards(RA_ART, DEC_ART, xy[0], xy[1],
                                                                                   rot=rot, parity=1.0)
            write_fits(os.path.join(d, name), cards)


def test_acquisition_frames_before_the_telescope_centres_are_no_drift(tmp_path):
    """Io 2017-12-13: the first frames 4.3' off, then centred all night. Measured from the first frames, the whole
    night looked displaced."""
    track_night(str(tmp_path), [(1675.0, 1372.0)] * 6 + [(1022.0, 1021.0)] * 200)
    assert tf.episodes(tf.night(str(tmp_path))) == []


def test_settling_then_steady_is_no_drift(tmp_path):
    """Io 2018-11-16: centred for six frames, moved 1.9' over four, then steady there for the rest of the block."""
    track_night(str(tmp_path), [(1010.0, 1032.0)] * 6 + [(1010 + 80 * k, 1032 - 75 * k) for k in range(1, 5)] +
                [(1339.0, 734.0)] * 150)
    assert tf.episodes(tf.night(str(tmp_path))) == []


def test_guider_bringing_the_field_back_after_a_pause_is_no_drift(tmp_path):
    """Io 2024-12-19 Sp0642-6707: after a 4.5-minute pause the field came back 1.3' off, and the guider walked it home
    in 40-px steps within six minutes. A jump without a pause, and one that stays, are still drifts."""
    home = [(1202.0, 1154.0)] * 60
    back = [(1429.0, 1173.0)] * 2 + [(1385.0, 1209.0)] * 3 + [(1350.0, 1165.0)] * 3 + [(1312.0, 1150.0)] * 3 + \
        [(1272.0, 1149.0)] * 3 + [(1232.0, 1149.0)] * 3
    d = str(tmp_path / 'a')
    os.makedirs(d)
    track_night(d, home + [None] * 8 + back + home)     # None: no frame solved... still frames, so no pause
    assert len(tf.episodes(tf.night(d))) == 1
    d = str(tmp_path / 'b')
    os.makedirs(d)
    t = dt.datetime(2024, 12, 19, 4, 0, 0)
    for i, xy in enumerate(home + back + home):
        ti = t + dt.timedelta(seconds=20 * i + (270 if i >= 60 else 0))       # the pause: 4.5 min after frame 60
        name = maxim_frame(d, i, ti)
        write_fits(os.path.join(d, name), list(tf.read_header(os.path.join(d, name)).items()) +
                   wcs_cards(RA_ART, DEC_ART, xy[0], xy[1], rot=0.6, parity=1.0))
    assert tf.episodes(tf.night(d)) == []


def test_scattered_wrong_solutions_are_unknown_not_drift(tmp_path):
    """Callisto 2024: solutions of frames 17 s apart land hundreds of pixels apart, inside the orientation and
    scale filters. The solutions are moving, not the pointing: the block cannot be judged."""
    import random
    rnd = random.Random(1)
    track_night(str(tmp_path), [(980 + rnd.uniform(-450, 450), 910 + rnd.uniform(-450, 450)) for _ in range(120)])
    r = tf.night(str(tmp_path))['Sp1356+4342']
    assert r['episodes'] == [] and r['solved'] == 0 and r['noisy'] == 120
    assert tf.flags({'Sp1356+4342': r}) == ('unknown', 'unknown')


def test_a_run_of_scattered_solutions_inside_a_steady_block_is_no_drift(tmp_path):
    import random
    rnd = random.Random(2)
    bad = [(1023 + rnd.choice((-1, 1)) * rnd.uniform(300, 700), 1021 + rnd.uniform(-700, 700)) for _ in range(8)]
    track_night(str(tmp_path), [(1023.0, 1021.0)] * 100 + bad + [(1023.0, 1021.0)] * 100)
    r = tf.night(str(tmp_path))['Sp1356+4342']
    assert r['episodes'] == [] and r['solved'] == 208 and r['noisy'] == 0


def test_excursions_off_the_field_and_back(tmp_path):
    """Io 2017-07-05 trappist-1: groups of frames 9-17' away (the target off the detector), then back home."""
    home = [(1023.0, 1021.0)] * 100
    away = [(1860.0, -330.0)] * 6 + [(668.0, -531.0)] * 3 + [(2128.0, -1729.0)] * 3
    track_night(str(tmp_path), home + away + home)
    (e,) = tf.night(str(tmp_path))['Sp1356+4342']['episodes']
    assert (e['i0'], e['i1'], e['n']) == (101, 112, 12)
    assert e['detail']['off_field_at'] == e['start']
    assert e['detail']['jump_arcmin'] == pytest.approx(math.hypot(837, 1351) * 0.35 / 60, abs=0.1)


def test_each_block_of_a_target_has_its_own_home(tmp_path):
    """A target observed twice in a night, another target between: pointed afresh, the second block's home can sit
    a few arcmin from the first's without anything drifting."""
    d = str(tmp_path)
    t = T0
    for i in range(30):
        astra_frame(d, i, t)
        t += dt.timedelta(minutes=1)
    for i in range(30):
        astra_frame(d, i, t, obj='Sp2135-3153', ra=323.81104, dec=-31.89607)
        t += dt.timedelta(minutes=1)
    for i in range(30):
        name = astra_frame(d, i, t, wcs=False)
        cards = [(k, v) for k, v in tf.read_header(os.path.join(d, name)).items()
                 if not k.startswith(('CTYPE', 'CRVAL', 'CRPIX'))]
        write_fits(os.path.join(d, name), cards + wcs_cards(206.10513, -40.3364, 752.0 + 600, 1038.0, rot=0.3))
        t += dt.timedelta(minutes=1)
    res = tf.night(d)
    assert tf.episodes(res) == [] and res['TOI-4336 c']['frames'] == 60 and res['TOI-4336 c']['solved'] == 60
    # without the night's other frames, a gap of BLOCK_GAP_MINUTES splits the blocks (31 min here)
    frames = sorted((f for f in tf.raw_frames(d) if f.object == 'TOI-4336 c'), key=lambda f: f.time)
    assert [len(b) for b in tf.target_blocks(frames)] == [30, 30]
    for f in frames[30:]:
        f.time -= dt.timedelta(minutes=20)
    assert [len(b) for b in tf.target_blocks(frames)] == [60]


def test_meridian_flip_is_a_new_reference_not_a_drift(tmp_path):
    """The field turns half a turn and the pointing lands elsewhere on the detector, all at once."""
    d = str(tmp_path)
    t = T0
    for i in range(40):
        flipped = i >= 20
        name = astra_frame(d, i, t, wcs=False)
        cards = tf.read_header(os.path.join(d, name))
        x, y = (2048 - 752.0, 2088 - 1038.0) if flipped else (752.0, 1038.0)
        write_fits(os.path.join(d, name), [(k, v) for k, v in cards.items() if not k.startswith(('CTYPE', 'CRVAL',
                                                                                                 'CRPIX'))] +
                   wcs_cards(206.10513, -40.3364, x, y, rot=180.3 if flipped else 0.3))
        t += dt.timedelta(minutes=2)
    assert tf.episodes(tf.night(d)) == []


def test_raw_header_and_log_line_solutions(tmp_path):
    """Astra frames carry their own WCS; frames without one take T7's CRPIX/CRVAL with a neighbour's matrix."""
    d = str(tmp_path)
    t, names = T0, []
    for i in range(30):
        names.append(astra_frame(d, i, t, wcs=i < 5))
        t += dt.timedelta(minutes=1)
    log = tmp_path / 'log.txt'
    lines = [t7_line(n, ok=True, crpix=(752.0 - (0 if i < 15 else 400.0), 1038.0), crval=(206.10513, -40.3364),
                     n=i + 1) for i, n in enumerate(names) if i >= 5]
    log.write_text('\n'.join(lines) + '\n')
    (e,) = tf.night(d, [str(log)])['TOI-4336 c']['episodes']
    assert e['kind'] == 'drift' and e['i0'] == 16 and e['detail']['jump_arcmin'] == pytest.approx(400 * 0.35 / 60,
                                                                                                    abs=0.05)


def test_wcs_cache_table(tmp_path):
    """The cache is astropy's ECSV: space-separated, double-quoted, the cards as one string of card images."""
    cards = ''.join(_card(k, v) for k, v in wcs_cards(209.172, 43.7164, 1052.0, 1017.0, parity=1.0)).rstrip()
    path = tmp_path / '20240211.ecsv'
    path.write_text('# %ECSV 1.0\n# ---\n# datatype:\n# - {name: raw, datatype: string}\n'
                    '# meta: !!omap\n# - {schema: 1}\n# schema: astropy-2.0\n'
                    'raw naxis1 naxis2 orchrel cards\n'
                    'Sp1356+4342-S001-R001-C001-I+z.fts 2046 2044 "" "' + cards + '"\n'
                    'bad-row 1 2\n')
    w = tf.wcs_cache(str(path))
    assert list(w) == ['Sp1356+4342-S001-R001-C001-I+z'] and w['Sp1356+4342-S001-R001-C001-I+z']['shape'] == (
        2044, 2046)
    x, y = tf.world_to_pixel(w['Sp1356+4342-S001-R001-C001-I+z'], 209.172, 43.7164)
    assert (x, y) == pytest.approx((1052.0, 1017.0))
    assert tf.wcs_cache(str(tmp_path / 'none.ecsv')) == {}
    assert tf.cache_path('/data', 'Artemis', '20240211') == '/data/PipelineOutput/wcs_cache/Artemis/20240211.ecsv'


def test_wcs_cache_written_by_astropy(tmp_path):
    table = pytest.importorskip('astropy.table')
    cards = ''.join(_card(k, v) for k, v in wcs_cards(209.172, 43.7164, 1052.0, 1017.0, parity=1.0)).rstrip()
    t = table.Table(rows=[['Sp1356+4342-S001-R001-C001-I+z.fts', 2046, 2044, cards, 'say "hi"']],
                    names=['raw', 'naxis1', 'naxis2', 'cards', 'proc'])
    t.meta['schema'] = 1
    t.write(str(tmp_path / 'c.ecsv'), format='ascii.ecsv')
    assert list(tf.wcs_cache(str(tmp_path / 'c.ecsv'))) == ['Sp1356+4342-S001-R001-C001-I+z']
