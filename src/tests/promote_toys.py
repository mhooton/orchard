"""
Toy data roots for the promotion step's tests and harness (tests/test_promote.py, tests/promote_harness.py).

A Toy writes what utils/promote.py and utils/promote_check.py read, at a size that takes milliseconds:
raw frame headers and schedule plan files under Observations/, and under PipelineOutput/v2 or v3 target
directories with aperture light curves promote_check can measure (a CATALOGUE with the target row first,
LIGHTCURVE_<ap> and ERROR_<ap> per star), *_output.fits, a 1/ run directory of frames with per-frame
sidecars and frame lists, lightcurves/ and stack files; night reductions with masters and image lists;
pipeline logs with the `date` lines and target lines ZLP_pipeline.sh writes; PDF reports.
"""
import datetime as dt
import hashlib
import os

import numpy as np
from astropy.io import fits

FILTER = 'I+z'


def gaia_id(name):
    """A fixed 19-digit stand-in Gaia ID for a toy star."""
    return '4' + str(int(hashlib.sha1(name.encode()).hexdigest(), 16))[:18]


def _header_only(path, **cards):
    hdr = fits.Header()
    for k, v in cards.items():
        hdr[k] = v
    fits.PrimaryHDU(header=hdr).writeto(path, overwrite=True)


def write_diff(path, ids, ap, frames=120, valid=None, noise=0.003, seed=0):
    """A *_diff.fits with one row per star (the first is the target) and `valid` finite target points."""
    rng = np.random.RandomState(seed)
    n = len(ids)
    flux = 1000.0 * (1.0 + noise * rng.standard_normal((n, frames)))
    if valid is not None:
        flux[0, valid:] = np.nan
    err = np.full((n, frames), 1000.0 * noise * 1.25)
    primary = fits.PrimaryHDU()
    primary.header['PIPE_V'] = 'v3'
    primary.header['DR3ID'] = ids[0]
    primary.header['DR2ID'] = 'None'
    cat = fits.BinTableHDU.from_columns([
        fits.Column(name='TARGET', format='J', array=np.array([1] + [0] * (n - 1))),
        fits.Column(name='GAIA_DR3_ID', format='25A', array=np.array(ids)),
        fits.Column(name='GAIA_DR2_ID', format='25A', array=np.array([''] * n))], name='CATALOGUE')
    fits.HDUList([primary, cat, fits.ImageHDU(flux, name='LIGHTCURVE_%d' % ap),
                  fits.ImageHDU(err, name='ERROR_%d' % ap)]).writeto(path, overwrite=True)


class Toy(object):
    def __init__(self, base):
        self.base = base
        self.po = os.path.join(base, 'PipelineOutput')

    # ---------------------------------------------------------------- Observations

    def raw(self, tel, night, target, n, filt=FILTER):
        d = os.path.join(self.base, 'Observations', tel, 'images', night)
        os.makedirs(d, exist_ok=True)
        stem = target.replace('--', ' ')
        for k in range(n):
            _header_only(os.path.join(d, '{}-{:04d}.fits'.format(target, k)), IMAGETYP='LIGHT', OBJECT=stem,
                         FILTER=filt)
        _header_only(os.path.join(d, 'bias-{}.fits'.format(target)), IMAGETYP='BIAS', OBJECT='bias', FILTER=filt)
        return d

    def plan(self, tel, night, target, gid):
        day = dt.datetime.strptime(night, '%Y%m%d').strftime('%Y-%m-%d')
        d = os.path.join(self.base, 'Observations', tel, 'schedule', 'Plans_by_date', day)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'Obj_{}.txt'.format(target.replace('--', ' '))), 'w') as f:
            f.write('; plan\n; toy\n; toy\n; {}\n#count 1\n{}\t01 00 00.00\t+10 00 00.00\n'.format(gid, target))

    # ---------------------------------------------------------------- PipelineOutput

    def target_dir(self, version, tel, night, target, gid=None, frames=120, valid=None, noise=0.003, aps=(3, 4, 5),
                   others=(), seed=0, n_frames_on_disk=6, filt=FILTER):
        """A target directory as T6-T10 leave it. others: further stars with their own light curves."""
        gid = gid or gaia_id(target)
        d = os.path.join(self.po, version, tel, 'output', night, target)
        run = os.path.join(d, '1')
        os.makedirs(os.path.join(d, 'lightcurves'), exist_ok=True)
        os.makedirs(run, exist_ok=True)
        stars = [gid] + list(others)
        for k, star in enumerate(stars):
            ids = [star] + [s for s in stars if s != star]
            for ap in aps:
                write_diff(os.path.join(d, '{}_{}_{}_{}_diff.fits'.format(star, filt, night, ap)), ids, ap,
                           frames=frames, valid=valid if k == 0 else None, noise=noise, seed=seed + 10 * k + ap)
            _header_only(os.path.join(d, '{}_{}_{}_output.fits'.format(star, filt, night)), NSTARS=len(stars))
            with open(os.path.join(d, 'lightcurves', '{}_{}_{}_5.png'.format(star, filt, night)), 'w') as f:
                f.write('png')
        _header_only(os.path.join(d, '{}_outstack_{}.fits'.format(target.upper(), filt)), STACK=1)
        _header_only(os.path.join(d, '{}_stack_catalogue_{}.fits'.format(target.upper(), filt)), CAT=1)
        names = []
        for k in range(n_frames_on_disk):
            stem = 'proc{}-{:04d}-{}'.format(target, k, filt)
            names.append(stem + '.fits')
            _header_only(os.path.join(run, stem + '.fits'), FRAME=k)
            for ext in ('.cat', '.fits.phot', '.fits.ell', '_psf.png'):
                with open(os.path.join(run, stem + ext), 'w') as f:
                    f.write(ext)
        for name in ('{}_processed.dat'.format(filt), '{}_stacked.dat'.format(filt)):
            with open(os.path.join(run, name), 'w') as f:
                f.write('\n'.join(names) + '\n')
        return d

    def reduction(self, version, tel, night, targets, masters=True, tag='new'):
        d = os.path.join(self.po, version, tel, 'output', night, 'reduction')
        os.makedirs(d, exist_ok=True)
        if masters:
            for m in ('1_MasterBias.fits', '1_MasterDark.fits', '1_MasterFlat_{}.fits'.format(FILTER)):
                _header_only(os.path.join(d, m), TAG=tag)
        for name in ('1_bias.list', '1_dark.list'):
            with open(os.path.join(d, name), 'w') as f:
                f.write(tag + '\n')
        for t in targets:
            with open(os.path.join(d, '1_image_{}.list'.format(t)), 'w') as f:
                f.write(tag + '\n')
        return d

    def log(self, version, tel, night, start, targets, specified=False, end=None, complete=True, stages=None):
        """A ZLP_pipeline.sh log: `date` lines in UTC, the 'Stages run:' line (given stages), the target line, the
        end."""
        d = os.path.join(self.po, version, tel, 'logs')
        os.makedirs(d, exist_ok=True)
        name = '{}_1_v3{}.log'.format(night, '_' + '_'.join(targets) if specified else '')
        end = end or start + dt.timedelta(hours=1)
        fmt = '%a %b %d %H:%M:%S UTC %Y'
        lines = ['Using python: /usr/local/bin/python', 'Environment set up', start.strftime(fmt)]
        if stages:
            lines.append('Stages run: ' + ' '.join('T{}'.format(n) for n in stages))
        lines += ['Using 20 cores', 'START T1',
                 'Using {} targets: {}'.format('specified' if specified else 'discovered', ' '.join(targets))]
        lines += ['TARGET: {}'.format(t) for t in targets]
        lines += ['START frame packing (int16)', end.strftime(fmt)]
        if complete:
            lines.append('PIPELINE COMPLETE')
        path = os.path.join(d, name)
        with open(path, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return path

    def pdf(self, version, tel, night):
        d = os.path.join(self.po, version, tel, 'reports')
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, '{}_{}.pdf'.format(tel, night))
        with open(path, 'w') as f:
            f.write('%PDF toy report {} {}\n'.format(tel, night))
        return path

    # ---------------------------------------------------------------- a whole night

    def night(self, tel, night, targets, start, version='v3', specified=None, reduction=True, pdf=True, **kw):
        """Raw frames, plan files and a run's output for each target, as one pipeline run leaves them.
        specified: the targets the run was given (a single-target run), else it discovered them all."""
        for t in targets:
            self.raw(tel, night, t, kw.get('frames', 120) + 1)
            self.plan(tel, night, t, kw.get('gids', {}).get(t) or gaia_id(t))
        ran = specified or targets
        self.log(version, tel, night, start, ran, specified=bool(specified))
        for t in ran:
            self.target_dir(version, tel, night, t, gid=kw.get('gids', {}).get(t), frames=kw.get('frames', 120),
                            seed=kw.get('seed', 0))
        if reduction:
            self.reduction(version, tel, night, targets)
        if pdf:
            self.pdf(version, tel, night)


def snapshot(root):
    """Every path under root with its size, mtime and inode: what a dry run must leave exactly as it was."""
    out = {}
    for dirpath, dirs, names in os.walk(root):
        for name in dirs + names:
            p = os.path.join(dirpath, name)
            st = os.lstat(p)
            out[os.path.relpath(p, root)] = (st.st_size if not os.path.isdir(p) else 0, st.st_mtime_ns, st.st_ino)
    return out


def tree_lines(root, top='PipelineOutput', max_files=6):
    """A compact listing for the harness: directories, and up to max_files names in each."""
    lines = []
    base = os.path.join(root, top)
    for dirpath, dirs, names in os.walk(base):
        dirs.sort()
        rel = os.path.relpath(dirpath, root)
        depth = rel.count(os.sep)
        shown = sorted(names)[:max_files]
        more = len(names) - len(shown)
        lines.append('{}{}/{}'.format('  ' * depth, os.path.basename(dirpath),
                                      '  [{}{}]'.format(', '.join(shown), ', +{} more'.format(more) if more else '')
                                      if names else ''))
    return lines
