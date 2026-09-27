"""
remove_duplicates was rewritten from an O(n^2) greedy scan to a grid-hashed
one because it was eating 44 s of the 60 s per-frame plate-solve budget on
star-rich ANDOR fields, which timed out whole Artemis nights.

The rewrite must be behaviour-preserving: the selected stars feed the WCS
fit, so a different selection is a different solution. These tests compare
the current implementation against the original loop, kept below verbatim as
a reference.

Run with pytest, or directly:  python src/tests/test_remove_duplicates.py
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from astrom.pointer_wcs import remove_duplicates

WIDTH, HEIGHT = 2046, 2044


def remove_duplicates_reference(detections, width, height):
    """The original O(n^2) implementation, kept as the behavioural oracle."""
    if len(detections) <= 1:
        return detections

    center = np.array([width / 2, height / 2])
    positions = detections[:, :2]
    scales = detections[:, 2]

    center_distances = np.linalg.norm(positions - center, axis=1)
    center_scores = 1.0 / (1.0 + center_distances / max(width, height))
    scale_scores = 1.0 / scales

    priority_scores = center_scores * scale_scores
    sorted_indices = np.argsort(priority_scores)[::-1]
    sorted_detections = detections[sorted_indices]

    unique_detections = []
    for detection in sorted_detections:
        is_duplicate = False
        pos = detection[:2]
        for existing in unique_detections:
            if np.linalg.norm(pos - existing[:2]) < 10.0:
                is_duplicate = True
                break
        if not is_duplicate:
            unique_detections.append(detection)

    return np.array(unique_detections)


def cases():
    rng = np.random.default_rng(20260927)

    # Scattered detections at three smoothing scales, as find_stars_multiscale
    # produces them.
    for n in (50, 500, 2000, 4000):
        yield (f'scattered n={n}',
               np.column_stack([rng.uniform(0, WIDTH, n),
                                rng.uniform(0, HEIGHT, n),
                                rng.choice([1.0, 2.0, 3.0], n),
                                rng.uniform(0, 1e5, n)]))

    # Tight clusters: the same star found at every scale, which is the case
    # the suppression radius exists for.
    for n in (2000, 6000):
        k = n // 4
        bx, by = rng.uniform(0, WIDTH, k), rng.uniform(0, HEIGHT, k)
        xs = np.concatenate([bx + rng.normal(0, 3, k) for _ in range(4)])
        ys = np.concatenate([by + rng.normal(0, 3, k) for _ in range(4)])
        yield (f'clustered n={len(xs)}',
               np.column_stack([xs, ys,
                                rng.choice([1.0, 2.0, 3.0], len(xs)),
                                rng.uniform(0, 1e5, len(xs))]))

    # Degenerate and boundary cases.
    yield 'empty', np.zeros((0, 4))
    yield 'single', np.array([[100.0, 100.0, 1.0, 5.0]])
    yield 'coincident', np.array([[100.0, 100.0, 1.0, 5.0],
                                 [100.0, 100.0, 2.0, 5.0]])
    # Either side of the 10 px suppression radius, including across a grid
    # cell boundary.
    yield 'radius edges', np.array([[100.0, 100.0, 1.0, 9.0],
                                    [109.9, 100.0, 1.0, 8.0],
                                    [110.1, 100.0, 1.0, 7.0],
                                    [100.0, 109.5, 1.0, 6.0]])
    # Negative and off-image coordinates must not break the cell indexing.
    yield 'off image', np.array([[-5.0, -5.0, 1.0, 4.0],
                                 [-1.0, -1.0, 2.0, 3.0],
                                 [WIDTH + 20.0, HEIGHT + 20.0, 1.0, 2.0]])


def main():
    ok = True
    print('%-18s %-7s %-7s %-7s %-10s %-10s %s'
          % ('case', 'n', 'ref', 'new', 'ref(s)', 'new(s)', 'identical'))
    for label, det in cases():
        t0 = time.perf_counter()
        expected = remove_duplicates_reference(det.copy(), WIDTH, HEIGHT)
        t_ref = time.perf_counter() - t0

        t0 = time.perf_counter()
        actual = remove_duplicates(det.copy(), WIDTH, HEIGHT)
        t_new = time.perf_counter() - t0

        same = expected.shape == actual.shape and np.array_equal(expected, actual)
        ok &= same
        print('%-18s %-7d %-7d %-7d %-10.4f %-10.4f %s'
              % (label, len(det), len(expected), len(actual), t_ref, t_new,
                 'YES' if same else 'NO <<<'))

    # The speed-up is the point of the change, so assert it holds.
    rng = np.random.default_rng(1)
    n = 4000
    big = np.column_stack([rng.uniform(0, WIDTH, n), rng.uniform(0, HEIGHT, n),
                           rng.choice([1.0, 2.0, 3.0], n), rng.uniform(0, 1e5, n)])
    t0 = time.perf_counter(); remove_duplicates(big.copy(), WIDTH, HEIGHT)
    t_new = time.perf_counter() - t0
    fast = t_new < 1.0
    ok &= fast
    print()
    print(('PASS  ' if fast else 'FAIL  ')
          + f'{n} detections deduplicated in {t_new:.3f}s (budget 1.0s)')

    print()
    print('ALL PASS' if ok else 'FAILURES PRESENT')
    return 0 if ok else 1


def test_remove_duplicates():
    assert main() == 0


if __name__ == '__main__':
    sys.exit(main())
