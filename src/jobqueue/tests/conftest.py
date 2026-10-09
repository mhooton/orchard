import datetime as dt
import os
import sys

import pytest

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if SRC not in sys.path:
    sys.path.insert(0, SRC)

from jobqueue.config import load_config  # noqa: E402
from jobqueue.resources import Snapshot  # noqa: E402
from jobqueue.store import Store  # noqa: E402


class Clock:
    def __init__(self, t):
        self.t = t

    def __call__(self):
        return self.t

    def advance(self, minutes=0, seconds=0):
        self.t += dt.timedelta(minutes=minutes, seconds=seconds)


class FixedResources:
    def __init__(self, load5=1.0, mem=400.0, util=None):
        self.load5, self.mem = load5, mem
        self.util = {'sda': 5.0, 'sdb': 5.0} if util is None else util

    def snapshot(self):
        return Snapshot(self.load5, self.mem, dict(self.util), 128)


@pytest.fixture
def clock():
    return Clock(dt.datetime(2026, 10, 9, 14, 0, 0))


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    monkeypatch.delenv('ORCHARD_QUEUE_CONFIG', raising=False)
    monkeypatch.delenv('ORCHARD_QUEUE_ROOT', raising=False)
    monkeypatch.setenv('PYTHONPATH', SRC + os.pathsep + os.environ.get('PYTHONPATH', ''))
    base = tmp_path / 'data'
    src = tmp_path / 'src'
    base.mkdir()
    src.mkdir()
    return load_config(overrides={'queue_root': str(tmp_path / 'queue'), 'basedir': str(base), 'src_dir': str(src),
                                  'python': sys.executable, 'eso_env_file': str(tmp_path / 'none.env')})


@pytest.fixture
def store(cfg, clock):
    s = Store(os.path.join(cfg['queue_root'], 'queue.sqlite'), now_fn=clock)
    yield s
    s.close()
