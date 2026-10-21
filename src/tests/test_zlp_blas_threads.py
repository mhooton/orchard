"""
ZLP_pipeline.sh runs every stage with one BLAS/OpenMP thread per process.

orchard-server's systemd scope allows 4,096 tasks (threads) for the whole container, and OpenBLAS starts a pool of up to
64 threads in every process that uses numpy. On 2026-10-20/21 concurrent backlog jobs ran the container out: OpenBLAS
raised SIGINT when pthread_create failed, so T1 died at numpy import (exit 130), T7 pool workers died mid-chunk and left
Pool.map waiting until the queue's timeout (13 jobs), and createlists' ThreadPool could not start (jobs 1075, 1076).

Run:  pytest src/tests/test_zlp_blas_threads.py -q
"""
import os
import re
import subprocess

SRC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZLP = os.path.join(SRC, 'main', 'ZLP_pipeline.sh')
VARS = ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')


def zlp_lines():
    with open(ZLP) as f:
        return f.read().splitlines()


def test_nothing_runs_before_the_thread_counts_are_set():
    lines = zlp_lines()
    first = min(i for i, line in enumerate(lines) if line.startswith('export ') and any(v in line for v in VARS))
    for line in lines[:first]:
        s = line.strip()
        assert not s or s.startswith('#') or re.match(r'set -o \w+$', s), \
            'ZLP_pipeline.sh runs {!r} before setting the BLAS thread counts'.format(s)
    exported = {re.match(r'export (\w+)=', line).group(1) for line in lines[first:first + len(VARS)]}
    assert exported == set(VARS)


def run_exports(env):
    exports = '\n'.join(line for line in zlp_lines() if line.startswith('export ') and any(v in line for v in VARS))
    script = 'set -o nounset\n' + exports + '\necho ' + ' '.join('${%s}' % v for v in VARS) + '\n'
    clean = {k: v for k, v in os.environ.items() if k not in VARS}
    clean.update(env)
    return subprocess.run(['bash', '-c', script], env=clean, capture_output=True, text=True, check=True).stdout.split()


def test_one_thread_unless_the_caller_sets_a_count():
    assert run_exports({}) == ['1', '1', '1']
    assert run_exports({'OMP_NUM_THREADS': '4'}) == ['1', '4', '1']
