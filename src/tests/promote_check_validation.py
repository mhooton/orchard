#!/usr/bin/env python3
"""
promote_check_validation.py
===========================

Validation cases for utils/promote_check.py, built from the 2026-10-09 audit
of the T12 sweep (mh_scratch/sweep_audit: v2_old_targets_classified.csv,
lc_compare.csv, v3_targets_classified.csv), and a score of its verdicts.

Historical sweeps: every v2_old target layer that has a v2 replacement.  The
candidate is the light curve the sweep put in v2, the incumbent the one it
moved to v2_old (layer N lives at v2_old/.../<target> nested N times).  As in
the audit, each old layer is compared with the current v2 directory.
Labels, using the audit's own definitions:
    regression: ap5 lost       old had an aperture-5 light curve, v2 has none (4)
    regression: all-NaN        v2's target light curve has no finite point (2)
    regression: points lost    v2 has under 90% of the old finite points (2)
    regression: noisier        v2's p2p over 1.2x the old, old not degenerate (5)
    improved: old degenerate   old p2p under 1e-3, the March 2026 Callisto reruns
    improved: identity         the sweep changed the Gaia ID (Wolf 359)
    improved: old had no ap5   old had no aperture-5 light curve
    fine: better / fine: similar   the rest of lc_compare, split at p2p ratio 1/1.2

Stranded: v3 target directories that pass T12's gate, or have aperture 5
without aperture 4, compared with whatever v2 holds at the same path (93).

Usage (inside the container):
    python tests/promote_check_validation.py pairs <sweep_audit dir> pairs.csv
    python utils/promote_check.py pairs pairs.csv --out verdicts.csv --raw-cache raw.json
    python tests/promote_check_validation.py score verdicts.csv
"""

import csv
import os
import sys
from collections import Counter, defaultdict

PO = '/data/SPECULOOSPipeline/PipelineOutput'
EXPECT = {'regression': ('HOLD', 'REJECT'), 'stranded': ('PROMOTE',),
          'improved': ('PROMOTE',), 'fine': ('PROMOTE',)}
SEVERITY = ['regression: ap5 lost', 'regression: all-NaN', 'regression: points lost', 'regression: noisier',
            'improved: old degenerate', 'fine: better', 'fine: similar']


def _f(x):
    try:
        return float(x)
    except ValueError:
        return float('nan')


def lc_label(r):
    old_n, new_n, old_p, new_p = _f(r['old_n']), _f(r['new_n']), _f(r['old_p2p']), _f(r['new_p2p'])
    if old_n > 0 and new_n == 0:
        return 'regression: all-NaN'
    if old_n > 0 and new_n < 0.9 * old_n:
        return 'regression: points lost'
    if old_p < 1e-3:
        return 'improved: old degenerate'
    if new_p > 1.2 * old_p:
        return 'regression: noisier'
    if new_p < old_p / 1.2:
        return 'fine: better'
    return 'fine: similar'


def build_pairs(audit, out):
    lc = defaultdict(list)
    with open(os.path.join(audit, 'lc_compare.csv'), newline='') as f:
        for r in csv.DictReader(f):
            lc[(r['tel'], r['datedir'], r['entry'], r['layer'])].append(lc_label(r))
    rows = []
    with open(os.path.join(audit, 'v2_old_targets_classified.csv'), newline='') as f:
        for r in csv.DictReader(f):
            if r['v2_exists'] != 'True' or r['testdir'] == 'True':
                continue
            key = (r['tel'], r['datedir'], r['entry'], r['layer'])
            if r['cls'].startswith('D'):
                label = 'regression: ap5 lost'
            elif r['cls'].startswith('C'):
                label = 'improved: identity'
            elif r['cls'].startswith('A'):
                label = 'improved: old had no ap5'
            elif lc[key]:
                label = min(lc[key], key=SEVERITY.index)
            else:
                label = 'other: %s' % r['cls']
            old = os.path.join(PO, 'v2_old', r['tel'], 'output', r['datedir'], *([r['entry']] * (1 + int(r['layer']))))
            new = os.path.join(PO, 'v2', r['tel'], 'output', r['datedir'], r['entry'])
            rows.append([r['tel'], r['datedir'], r['entry'], new, old, label])
    with open(os.path.join(audit, 'v3_targets_classified.csv'), newline='') as f:
        names = {'1': 'stranded: never published', '2': 'stranded: newer than v2', '4': 'stranded: ap5 without ap4'}
        for r in csv.DictReader(f):
            c = r['cls'][:1]
            if c in names and r['testdir'] != 'True' and r['layer'] == '0':
                rows.append([r['tel'], r['datedir'], r['entry'],
                             os.path.join(PO, 'v3', r['tel'], 'output', r['datedir'], r['entry']),
                             os.path.join(PO, 'v2', r['tel'], 'output', r['datedir'], r['entry']), names[c]])
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['tel', 'date', 'target', 'candidate', 'incumbent', 'label'])
        w.writerows(rows)
    print(Counter(r[-1] for r in rows))


def score(path):
    with open(path, newline='') as f:
        rows = list(csv.DictReader(f))
    print('%-32s %8s %8s %8s   expected' % ('label', 'PROMOTE', 'HOLD', 'REJECT'))
    for lab in sorted({r['label'] for r in rows}):
        c = Counter(r['verdict'] for r in rows if r['label'] == lab)
        print('%-32s %8d %8d %8d   %s' % (lab, c['PROMOTE'], c['HOLD'], c['REJECT'],
                                          '/'.join(EXPECT.get(lab.split(':')[0], ('?',)))))
    wrong = [r for r in rows if r['verdict'] not in EXPECT.get(r['label'].split(':')[0], ())]
    print('\n%d not as expected:' % len(wrong))
    for r in sorted(wrong, key=lambda r: (r['label'], r['tel'], r['date'])):
        print('  %-26s %-8s %s %-16s %-7s %s' % (r['label'], r['tel'], r['date'], r['target'], r['verdict'],
                                                 r['reasons'] or r['notes']))


if __name__ == '__main__':
    if sys.argv[1] == 'pairs':
        build_pairs(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == 'score':
        score(sys.argv[2])
