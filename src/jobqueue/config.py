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

    # code releases (releases.py): jobs run from <releases_dir>/current, a frozen copy of the deployed tree
    'releases_dir': None,                     # default <queue_root>/releases
    # left out of a release: caches, logs and core dumps that live in the production tree; '/x' anchors at the root
    'release_exclude': ['__pycache__', '.pytest_cache', '.git', '.idea', 'catcache', 'core', 'core.[0-9]*', '._*',
                        '.DS_Store', '*.pyc', '*.log', '*.release-tmp', '/1'],
    # not copied into a release but symlinked to the live tree's file: credentials (reporting/.env), so each keeps
    # its owner-only mode, exists once, and no sync or rollback writes an old one back. Patterns as release_exclude
    'release_link': ['.env'],
    'test_pylib': '/data/SPECULOOSPipeline/mh_scratch/k1-deploy-check/pylib',   # pytest for the container
    'test_paths': ['tests', 'jobqueue/tests'],

    # dispatcher
    'poll_seconds': 45,
    'max_cores': 112,
    'load_max': 100.0,
    'mem_min_gb': 64.0,
    'disk_util_max': 70.0,
    'disk_devices': ['sda', 'sdb'],           # sda = PipelineOutput, sdb = data/raw on appct
    'max_disk_heavy_starts_per_poll': 1,      # %util needs a poll to show a new job's I/O
    # cores kept free for P0 in the nightly window; 0: a backlog job at the head of the queue runs whenever the
    # capacity is there, and a P0 job that finds it full waits for the next job to finish
    'p0_reserve_cores': 0,
    'p0_reserve_window_utc': ['11:00', '01:00'],
    # one pipeline job per telescope at a time. Off since each wcsfit call got its own catcache and the report
    # tables a lock (K7); set it to go back to one at a time (pipeline jobs and external ZLP runs keep tel: locks)
    'telescope_lock': False,
    # concurrent ESO download/fetch jobs; P2/P3 fetches (backfill, q add --download) also take eso_bulk, so they
    # never hold both eso slots and the nightly downloads always find one free
    'semaphores': {'eso': 2, 'eso_bulk': 1},
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

    # the add-only fetch (fetch.py): downloads in flight per job, tries per file, first backoff (doubles each
    # try), consecutive network failures before a job gives up for now, rows per staging batch, and the
    # download rate assumed for estimates (MB/s of ESO's compressed size)
    'fetch': {'workers': 3, 'attempts': 4, 'backoff_seconds': 30, 'breaker': 8, 'batch': 25, 'mb_per_s': 4.0},

    # the promotion step (utils/promote.py) that replaces T12's sweep: a 'promote' job after a pipeline job run with
    # --no_T12 (q add --promote) judges the night's v3 output with promote_check and publishes what passes into v2
    'promote': {
        'db': None,                    # the journal of verdicts and promotions; default <queue_root>/promote.sqlite
        'raw_cache': None,             # raw-frame counts per night (check A4); default <queue_root>/promote_raw.sqlite
        'identity': 'plan',            # A3: 'plan' (schedule plan files), or 'table': <basedir>/ml_40pc_v2.csv for Sp
        'thresholds': {},              # overrides of promote_check.THRESHOLDS, e.g. {"max_scatter_ratio": 1.3}
        'cores': 2,
        'est_minutes': 15.0,
        'priority': 1000000000.0,      # ahead of the pipeline jobs of its class: it is short and holds up their results
        # also run the watcher's, the look-back's and fetch reruns' pipeline jobs with --no_T12 and a promote job
        # after each; off until the nightly runs move from T12 to promotion
        'nightly': False,
    },

    'lookback': {
        'nights': 30,
        'max_fetch_rows_per_night': 3000,
        'max_fetch_rows_per_run': 10000,
        'max_reruns_per_night': 2,
        'max_fetches_per_night': 3,
    },

    # the notifier (notify.py, digest.py): the morning digest and incident alerts, the pipeline's only email.
    # 'outbox' writes each message as an .eml file under <queue_root>/notify/outbox and sends nothing; 'smtp'
    # also sends it, to `to`, with EMAIL_ADDRESS, EMAIL_PASSWORD and SMTP_SERVER read from env_file
    'notify': {
        'transport': 'outbox',
        'env_file': None,
        'to': [],
        'subject_prefix': '[orchard]',
        'heartbeat_url': None,           # fetched after every digest: an outside check that alerts if it stops
        'late_after_utc': '21:00',       # transfer-log frames still missing at ESO at this time the next day
        'dark_min_nights': 3,            # nights a telescope sits dark while two others observe before an incident
        'dark_lookback_nights': 120,
        'dark_min_frames': 150,          # fewer frames on disk and no light curve: calibrations only, not observed
        'chilean_min_fraction': 0.5,     # share of the planned night in SPOCK's ch_ blocks that makes it Chilean
        'failed_lookback_nights': 3,     # P0 failures on nights this recent open incidents
        'signature_alert_count': 10,     # failures with one signature in 24 h that open an incident
        'stall_hours': 6.0,              # jobs due, none running and none started for this long
        'dispatcher_stale_minutes': 30,
        'disk_min_free_tb': 2.0,
        'disk_min_free_pct': 5.0,
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
        'quarantine': os.path.join(root, 'quarantine'),
        'manifest': os.path.join(root, 'fetch', 'manifest.csv'),
        'reports': os.path.join(root, 'reports'),
        'outbox': os.path.join(root, 'notify', 'outbox'),
        'promote_db': (cfg.get('promote') or {}).get('db') or os.path.join(root, 'promote.sqlite'),
        'promote_raw': (cfg.get('promote') or {}).get('raw_cache') or os.path.join(root, 'promote_raw.sqlite'),
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
