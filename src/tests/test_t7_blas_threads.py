"""
T7 runs with one BLAS thread per process (astrom/astrometry.py, top of the module).

T7's parallelism is its process pool. On appct (2026-10-20/21) forked workers whose OpenBLAS had started its 64-thread
pool deadlocked: all 64 threads in futex_wait for hours, out of reach of the 60 s SIGALRM, until the queue's job
timeout. numpy and scipy each load their own OpenBLAS; both read OPENBLAS_NUM_THREADS / OMP_NUM_THREADS when they are
first loaded, so astrometry.py must set them before anything imports numpy.

The thread count needs Linux /proc: run inside orchard-server.
Run:  pytest src/tests/test_t7_blas_threads.py -q
"""
import ast
import os
import subprocess
import sys

import pytest

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASTROMETRY = os.path.join(SRC, 'astrom', 'astrometry.py')
VARS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')

# Imports the module as T7's process does, then does BLAS work big enough for OpenBLAS to use its pool, in numpy's
# library and in scipy's, and prints the process's thread count.
CHECK = r'''
import os, sys
sys.path.insert(0, sys.argv[1])
if sys.argv[2] == 'astrometry':
    import astrom.astrometry
import numpy as np
a = np.random.default_rng(1).random((600, 600))
np.dot(a, a)
try:
    import scipy.linalg
    scipy.linalg.inv(a)
except ImportError:
    pass
print(len(os.listdir('/proc/self/task')))
'''


def _threads(first_import, blas_threads):
    """Threads of a fresh process after BLAS work, with the thread variables set to blas_threads in its environment.
    Never OpenBLAS's default: under load on appct its 64 threads per process fail to start."""
    env = dict(os.environ)
    env.update({k: str(blas_threads) for k in VARS})
    out = subprocess.run([sys.executable, '-c', CHECK, SRC, first_import], env=env, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True, check=True, timeout=300)
    return int(out.stdout.split()[-1])


@pytest.mark.skipif(not os.path.isdir('/proc/self/task'), reason='counts threads in /proc: run inside orchard-server')
def test_t7_runs_with_one_blas_thread():
    pytest.importorskip('twirl')                     # astrometry.py imports it
    with_t7 = _threads('astrometry', blas_threads=4)         # astrometry.py overrides what the environment says
    assert with_t7 <= 3, '{} threads after importing astrom.astrometry'.format(with_t7)
    if (os.cpu_count() or 1) >= 4:                  # the same work without it starts the BLAS pools
        assert _threads('numpy', blas_threads=4) > with_t7


def test_nothing_imports_numpy_before_the_threads_are_set():
    """The assignment must come before every import but `import os`; an editor's import sorting would undo it."""
    tree = ast.parse(open(ASTROMETRY).read())
    for node in tree.body:
        if isinstance(node, ast.For) and 'OPENBLAS_NUM_THREADS' in ast.dump(node):
            break
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = [a.name for a in node.names] if isinstance(node, ast.Import) else [node.module]
            assert names == ['os'], 'imports {} before setting the BLAS threads'.format(names)
    else:
        pytest.fail('astrometry.py no longer sets the BLAS threads')
