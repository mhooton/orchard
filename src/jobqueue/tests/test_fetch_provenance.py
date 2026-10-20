"""The queue's fetch writes raw frames as docs/versioning.md ("Raw frames") says: a datacube's frames and Astra frames
are stamped at ingest and their checksums verify; any other frame is placed byte for byte as ESO delivered it."""
import io
import os

import numpy as np
from astropy.io import fits

from jobqueue import fetch, jobs
from utils import orchard_version as ov

ACP = 'SPECU1.2019-04-15T23:00:00.000'
ANDOR = 'SPECU2.2026-10-16T01:00:00.000'
CUBE = 'SPECU4.2026-10-16T02:00:00.000'


def fits_bytes(hdus):
    buf = io.BytesIO()
    fits.HDUList(hdus).writeto(buf, checksum=True)       # ESO's files carry CHECKSUM/DATASUM
    return buf.getvalue()


def frame(observer, shape=(6, 9), imagetyp='Light Frame'):
    h = fits.Header()
    h['IMAGETYP'], h['OBSERVER'], h['DATE-OBS'] = imagetyp, observer, '2026-10-16T01:00:00.000'
    data = np.random.default_rng(0).integers(200, 60000, size=shape, dtype=np.uint16)
    return fits_bytes([fits.PrimaryHDU(data=data, header=h)])


def cube(n=2):
    h = fits.Header()
    h['IMAGETYP'], h['OBSERVER'] = 'Light Frame', 'Astra'
    data = np.random.default_rng(1).integers(200, 16000, size=(n, 4, 1024), dtype=np.uint16)
    meta = fits.BinTableHDU.from_columns([
        fits.Column(name='FILENAME', format='40A', array=np.array(['SPECULOOS4.2026-10-16T02_00_{:02d}.000'.format(i)
                                                                   for i in range(n)])),
        fits.Column(name='JD-OBS', format='D', array=2461329.5 + np.arange(n) / 86400.0),
    ], name='METADATA')
    return fits_bytes([fits.PrimaryHDU(data=data, header=h), meta])


def fake_fetcher(blobs):
    class Fetcher:
        """Stands in for EsoFileFetcher: 'downloads' the blobs into the staging directory."""
        def __init__(self, token, **kwargs):
            self.failed, self.retry_later, self.bytes = {}, set(), 0

        def __call__(self, dp_ids, dest):
            out = []
            for d in dp_ids:
                with open(os.path.join(dest, d + '.fits'), 'wb') as f:
                    f.write(blobs[d])
                out.append(os.path.join(dest, d + '.fits'))
            return out
    return Fetcher


def verifies(path):
    with fits.open(path, do_not_scale_image_data=True) as h:
        return h[0].verify_datasum() == 1 and h[0].verify_checksum() == 1


def test_fetch_stamps_what_it_writes_and_leaves_downloads_as_delivered(cfg, monkeypatch, tmp_path):
    blobs = {ACP: frame('SPECULOOS'), ANDOR: frame('Astra'), CUBE: cube()}
    monkeypatch.setattr(fetch, 'EsoFileFetcher', fake_fetcher(blobs))
    monkeypatch.setattr(jobs, 'eso_token_fn', lambda archive: None)
    dest = tmp_path / 'trial'
    res = jobs.fetch_rows(cfg, None, 'Callisto', '20261015', [{'dp_id': d} for d in blobs], dest=str(dest),
                          log=lambda m: None)
    assert res['added'] == 4 and not res['invalid']
    v = ov.version()
    # downloaded and kept as ESO delivered it: byte for byte, no stamp
    assert open(dest / (ACP + '.fits'), 'rb').read() == blobs[ACP]
    assert 'ORCHVER' not in fits.getheader(dest / (ACP + '.fits'))
    # rotated at ingest
    h = fits.getheader(dest / (ANDOR + '.fits'))
    assert h['ORCHVER'] == v and h['ASTRAROT'] == 90 and 'ORCHMOD' not in h
    assert [str(c) for c in h['HISTORY'] if str(c).startswith('orchard ')] == ['orchard {} astra_transform'.format(v)]
    assert verifies(dest / (ANDOR + '.fits'))
    # unpacked, then mirrored
    for i in range(2):
        path = dest / 'SPECU4.2026-10-16T02:00:{:02d}.000.fits'.format(i)
        h = fits.getheader(path)
        assert h['ORCHVER'] == v and h['ASTRAMIR'] == 'Y'
        assert [str(c) for c in h['HISTORY'] if str(c).startswith('orchard ')] == [
            'orchard {} unpack_datacubes'.format(v), 'orchard {} astra_transform'.format(v)]
        assert verifies(path)
