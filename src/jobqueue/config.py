"""Queue configuration: built-in defaults, optionally overridden by a JSON file.

The file is found from --config, then $ORCHARD_QUEUE_CONFIG, then <queue_root>/config.json.
$ORCHARD_QUEUE_ROOT overrides queue_root. Only keys present in the file change; nested
dictionaries are merged.
"""
import copy
import json
import os

DEFAULTS = {
    # 'live' launches jobs; 'shadow' only records what it would launch and simulates the run.
    'mode': 'live',
    'queue_root': '/data/SPECULOOSPipeline/queue',
    'basedir': '/data/SPECULOOSPipeline',
    'src_dir': '/opt/orchard/src',
    'eso_env_file': '/opt/orchard/src/reporting/.env',
    'python': 'python',

    # dispatcher
    'poll_seconds': 45,
    'max_cores': 112,
    'load_max': 100.0,
    'mem_min_gb': 64.0,
    'disk_util_max': 70.0,
    'disk_devices': ['sda', 'sdb'],           # sda = PipelineOutput, sdb = data/raw on appct
    'max_disk_heavy_starts_per_poll': 1,      # %util needs a poll to show a new job's I/O
    'p0_reserve_cores': 40,
    'p0_reserve_window_utc': ['11:00', '01:00'],
    'telescope_lock': True,                   # until a per-run catcache lands (decision K7)
    'semaphores': {'eso': 2},                 # concurrent ESO download/fetch jobs
    'timeout_factor': 3.0,
    'min_timeout_minutes': 120,
    'kill_grace_seconds': 120,
    'retry_delay_minutes': 60,
    'max_retries': 2,
    'vanished_requeues': 1,
    'heartbeat_stale_seconds': 600,
    'transient_exit_codes': [75, 137, 143, -9, -15],
    'transient_patterns': [
        r'No space left on device',
        r'Connection (reset|refused|aborted)',
        r'Temporary failure in name resolution',
        r'Name or service not known',
        r'(Read|Connect|Connection) ?timed out|ReadTimeout|ConnectTimeout',
        r'Max retries exceeded',
        r'RemoteDisconnected|Remote end closed',
        r'Stale file handle',
        r'Input/output error',
        r'Cannot allocate memory|MemoryError',
        r'Resource temporarily unavailable',
        r'Too many open files',
        r'database is locked',
        r'Broken pipe',
        r'50[234] (Server Error|Service Unavailable|Bad Gateway|Gateway Time-out)',
    ],

    # what a nightly pipeline job looks like
    'pipeline_cores': 20,
    'download_cores': 2,
    'pipeline_runname': '1',
    'pipeline_cthresh': '8',
    'pipeline_sthresh': '2',
    'minutes_per_1000_frames': {'ANDOR': 69.0, 'SPIRIT': 22.0},
    'min_estimate_minutes': 15.0,
    # Callisto carried SPIRIT on these nights (measured; monitoring.py's dates are wrong)
    'spirit_spells': [['20220111', '20220111'], ['20220115', '20220115'],
                      ['20220510', '20230315'], ['20250226', None]],

    'telescopes': {
        'Io': {'source': 'eso', 'prog_id': '60.A-9009(A)', 'force_platesolve': True},
        'Europa': {'source': 'eso', 'prog_id': '60.A-9009(B)', 'force_platesolve': True},
        'Ganymede': {'source': 'eso', 'prog_id': '60.A-9009(C)', 'force_platesolve': True},
        'Callisto': {'source': 'eso', 'prog_id': '60.A-9009(D)', 'force_platesolve': True},
        'Artemis': {'source': 'raw_dir', 'force_platesolve': False},
    },

    'watcher': {
        'watch_days': 2,               # last night and the one before
        'stable_minutes': 25,          # same count on two polls at least this far apart
        'quiet_minutes': 20,           # and nothing modified at ESO for this long
        'partial_stable_minutes': 180, # below the transfer-log count: wait this long before giving up
        'no_data_note_after_utc': '21:00',
        'artemis_manifest': 'Data_Download.txt',
        'artemis_deadline_utc': '22:00',
        'recent_days': 5,              # SSO_download.py is only ever run on nights this recent
    },

    'lookback': {
        'nights': 30,
        'max_fetch_rows_per_night': 3000,
        'max_fetch_rows_per_run': 10000,
        'max_reruns_per_night': 2,
        'max_fetches_per_night': 3,
    },

    # Keep a lower-class job off a telescope whose nightly job is due before it would finish.
    'p0_guard': {
        'enabled': True,
        'windows_utc': {'Io': ['13:00', '17:00'], 'Europa': ['13:00', '17:00'],
                        'Ganymede': ['13:00', '17:00'], 'Callisto': ['13:00', '17:00'],
                        'Artemis': ['10:00', '13:00']},
    },
}


def _merge(base, extra):
    for k, v in extra.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _merge(base[k], v)
        else:
            base[k] = v
    return base


def load_config(path=None, overrides=None):
    cfg = copy.deepcopy(DEFAULTS)
    if os.getenv('ORCHARD_QUEUE_ROOT'):
        cfg['queue_root'] = os.environ['ORCHARD_QUEUE_ROOT']
    path = path or os.getenv('ORCHARD_QUEUE_CONFIG')
    if not path:
        candidate = os.path.join(cfg['queue_root'], 'config.json')
        path = candidate if os.path.exists(candidate) else None
    if path:
        with open(path) as f:
            _merge(cfg, json.load(f))
        cfg['_config_path'] = path
    if overrides:
        _merge(cfg, overrides)
    if cfg['mode'] not in ('live', 'shadow'):
        raise ValueError('mode must be live or shadow, not {!r}'.format(cfg['mode']))
    return cfg


def queue_paths(cfg):
    root = cfg['queue_root']
    return {
        'root': root,
        'db': os.path.join(root, 'queue.sqlite'),
        'lock': os.path.join(root, 'dispatcher.lock'),
        'logs': os.path.join(root, 'logs'),
        'job_logs': os.path.join(root, 'logs', 'jobs'),
        'fetch': os.path.join(root, 'fetch'),
        'staging': os.path.join(root, 'staging'),
        'reports': os.path.join(root, 'reports'),
    }


def ensure_dirs(cfg):
    paths = queue_paths(cfg)
    for key in ('root', 'logs', 'job_logs', 'fetch', 'staging', 'reports'):
        os.makedirs(paths[key], exist_ok=True)
    return paths


def eso_telescopes(cfg):
    return [t for t, v in cfg['telescopes'].items() if v.get('source') == 'eso']


def camera(cfg, telescope, night):
    """'SPIRIT' or 'ANDOR' for timing estimates."""
    if telescope != 'Callisto':
        return 'ANDOR'
    for a, b in cfg['spirit_spells']:
        if a <= night and (b is None or night <= b):
            return 'SPIRIT'
    return 'ANDOR'
