"""The LHATPRO PWV fetch (T10): where it goes, and what the log says when it cannot get there.

From 2026-10-18 every name lookup inside orchard-server intermittently failed (the host's first
nameserver answered SERVFAIL and docker 1.13's embedded DNS passes that on). getLHATPROdata raised
ConnectionError, add_water_vapour's bare except printed only "No PWV data for: <window>", and the
light curves went out without a PWV correction while the logs read like an empty archive window.

These tests pin the query to https with a timeout, make HTTP errors raise instead of parsing an
error page, and make the "No PWV data" line carry the reason.

Run with `pytest src/tests/test_pwv_lhatpro.py`.
"""

import importlib
import os
import re
import sys
import types

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

requests = pytest.importorskip('requests')

# What the wdb returns: tab-separated, every field space-padded, blanks as all spaces.
TSV = (
    "Platform\tDate time\tIntegration time [s]\tIR temperature [Celsius]\t"
    "Liquid water path [g/m**2]\tPrecipitable Water Vapour [mm]\n"
    "A                             \t2026-10-07T00:00:01 \t110                 \t-95.588                       \t"
    "                              \t                              \n"
    "C                             \t2026-10-07T00:02:14 \t103                 \t-81.032                       \t"
    "-11.147                       \t0.871                         \n"
    "C                             \t2026-10-07T00:04:14 \t103                 \t-80.990                       \t"
    "-11.120                       \t0.902                         \n"
)


def _pwvgrid():
    try:
        return importlib.import_module('water_vapour.pwvGrid')
    except ImportError as e:      # e.g. scipy >= 1.14 dropped scipy.integrate.simps
        pytest.skip('water_vapour.pwvGrid does not import here: {}'.format(e))


class _Response:
    def __init__(self, text='', status=200):
        self.text = text
        self.status_code = status

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError('{} Server Error'.format(self.status_code), response=self)


def test_query_goes_to_https_with_a_timeout(monkeypatch):
    pwv = _pwvgrid()
    calls = []

    def post(url, **kwargs):
        calls.append((url, kwargs))
        return _Response(TSV)

    monkeypatch.setattr(pwv.requests, 'post', post)
    df = pwv.getLHATPROdata('2026-10-07T00:00:00', '2026-10-07T01:00:00', format='jd')

    (url, kwargs), = calls
    assert url == 'https://archive.eso.org/wdb/wdb/asm/lhatpro_paranal/query'
    assert kwargs.get('timeout'), 'an unanswered POST would hang T10'
    assert kwargs['data']['start_date'] == '2026-10-07T00:00:00..2026-10-07T01:00:00'
    assert list(df['pwv']) == [0.871, 0.902]
    assert df.index[0] == pytest.approx(2461320.5 + (2 * 60 + 14) / 86400.0, abs=1e-6)


def test_http_error_raises_rather_than_parsing_the_error_page(monkeypatch):
    pwv = _pwvgrid()
    monkeypatch.setattr(pwv.requests, 'post', lambda url, **kw: _Response('<html>503</html>', 503))
    with pytest.raises(requests.HTTPError):
        pwv.getLHATPROdata('2026-10-07T00:00:00', '2026-10-07T01:00:00', format='jd')


def _splightcurve():
    try:
        return importlib.import_module('diff_photometry.SPlightcurve')
    except ImportError:
        pass
    saved = dict(sys.modules)
    try:
        for name in ('fitsio', 'water_vapour', 'water_vapour.water_vapour', 'water_vapour.pwvGrid'):
            try:
                importlib.import_module(name)
            except ImportError:
                stub = types.ModuleType(name)
                stub.getLHATPROdata = None
                stub.FITS = object
                sys.modules[name] = stub
        sys.modules.pop('diff_photometry.SPlightcurve', None)
        try:
            return importlib.import_module('diff_photometry.SPlightcurve')
        except ImportError as e:
            pytest.skip('SPlightcurve does not import here: {}'.format(e))
    finally:
        for name in list(sys.modules):
            if name not in saved and not name.startswith('diff_photometry'):
                sys.modules.pop(name, None)


def test_unreachable_archive_is_reported_with_its_reason(monkeypatch, capsys):
    splc = _splightcurve()
    calls = []

    def unreachable(dateS, dateE, **kw):
        calls.append((dateS, dateE))
        raise requests.ConnectionError("HTTPSConnectionPool(host='archive.eso.org', port=443): Max retries "
                                       "exceeded (Failed to resolve 'archive.eso.org')")

    monkeypatch.setattr(splc, 'getLHATPROdata', unreachable)
    tel = types.SimpleNamespace()
    # T10 catches this around pwv_correct and writes the light curve uncorrected, logging "PWV Failed: <e>",
    # so the message is what reaches the log. It used to be numpy's "All-NaN slice encountered".
    with pytest.raises(ValueError) as err:
        splc.Telescope.add_water_vapour(tel, [2461320.6, 2461320.7], False)

    line, = [l for l in capsys.readouterr().out.splitlines() if l.startswith('No PWV data for:')]
    assert 'ConnectionError' in line and "Failed to resolve 'archive.eso.org'" in line
    assert len(calls) == 2                      # one retry, as before
    assert str(err.value).startswith('no LHATPRO PWV data (ConnectionError: ')
    assert 'All-NaN' not in str(err.value)
    assert tel.pwv == [] and tel.pwv_t == []


def test_an_empty_window_says_so(monkeypatch):
    splc = _splightcurve()
    monkeypatch.setattr(splc, 'getLHATPROdata', lambda dateS, dateE, **kw: None)
    with pytest.raises(ValueError, match=r'no LHATPRO PWV data \(the archive returned no valid rows\)'):
        splc.Telescope.add_water_vapour(types.SimpleNamespace(), [2461320.6, 2461320.7], False)


def test_no_plain_http_to_the_eso_archive():
    """http://archive.eso.org is a redirect to https; nothing in the pipeline should rely on it."""
    offenders = []
    for root, dirs, files in os.walk(SRC):
        dirs[:] = [d for d in dirs if d not in ('tests', '__pycache__')]
        for fn in files:
            if fn.endswith(('.py', '.sh')):
                path = os.path.join(root, fn)
                with open(path, errors='replace') as fh:
                    for n, line in enumerate(fh, 1):
                        if re.search(r'http://archive\.eso\.org', line):
                            offenders.append('{}:{}'.format(os.path.relpath(path, SRC), n))
    assert not offenders, offenders


if __name__ == '__main__':
    sys.exit(pytest.main([__file__, '-v']))
