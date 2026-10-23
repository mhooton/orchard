# Operating rules for Claude sessions in orchard

orchard runs in production on appct, inside the `orchard-server` container, driven by the job queue
([docs/job-queue.md](docs/job-queue.md)). These rules exist because each one has been broken once. Write in
British English.

## Production

- **Deploy only with `q release deploy`, and only when Matthew explicitly asks for that deploy.** Approval of a
  PR is not approval to deploy it. Rollback is `q release activate <previous release id>` (plus `q stop` if
  `jobqueue/` changed). See "Deploying code: releases" in docs/job-queue.md.
- **Never rsync or copy a `src/` tree into production** (`/appct/data/speculoos/orchard/src`, container
  `/opt/orchard/src`). It reverts `jobqueue/`, bypasses releases and edits the code every new job runs.
- **Never delete or move data** under the data root (`PipelineOutput/`, `Observations/`, `queue/`,
  `mh_scratch/`) without Matthew's explicit OK for that specific action.
- **Never edit `queue/config.json` or the crontab** without an explicit instruction; back up first, and say
  where the backup is.
- **Never open `queue.sqlite` from appcs**: it is on NFS there. Use `q` on appct, or a read-only sqlite
  backup inside the container.
- **Launch once, check status separately.** A detached launch succeeds even when a later command in the same
  ssh string fails; a retry loop around it launches twice. Before relaunching anything, check that nothing is
  already running.

## Email

- **Never send email.** Drafts only, as files in `.claude/emails/` (or Outlook drafts that Matthew sends).
  The pipeline's notifier is the only sender, and it is configured by Matthew.
- **Tests must never reach production config or real transports.** Release deploys run the suite inside the
  production container, where `queue/config.json` exists. Fixtures set their own queue root; sending code
  takes fakes by default.

## Credentials

- **Never print any part of a `.env` file**, not even one line or its length. Check a key exists with
  `grep -c '^KEY=' file`. Never widen its permissions beyond 600.

## Running things on appct

- **Run orchard's Python only inside the container**: `docker exec orchard-server python ...` with container
  paths (`/data/SPECULOOSPipeline`, `/opt/orchard/src`, `/gaia_database`). Not appct's host python, not
  appcs's miniconda.
- **Paths:** appcs `/appct/data/SPECULOOSPipeline` = container `/data/SPECULOOSPipeline`. Stage scratch
  work under `mh_scratch/<job>/`.
- **Two ssh hops lose quotes.** From the Mac, `ssh appcs 'ssh appct docker exec orchard-server bash -c "cd X
  && python y.py"'` runs `python y.py` on the appct host. Pass only space-free arguments through the hops
  (`... orchard-server bash /abs/script.sh arg`), or put the commands in a script file.
- **`docker exec` needs `-i` to read a heredoc**; without it the command waits for ever. appct's docker is
  old: no `docker exec -w`; use absolute paths and `-e`.
- **`python -m pkg` runs from the working directory first**: the container's default cwd is `/opt/orchard/src`,
  so `python -m jobqueue` or `python -m pytest` runs the *deployed* package whatever `PYTHONPATH` says. Run from
  the snapshot directory (`bash -c "cd <snapshot> && python -m ..."`).
- **`pgrep -f` / `pkill -f` match their own shell** when the pattern is in the ssh command line. Bracket the
  first letter (`pgrep -af "[Z]LP_pipeline.sh"`), and never run `pkill -f` from a remote one-liner.
- **zsh does not word-split variables.** `H="ssh appcs"; $H cmd` fails on this Mac. Use a shell function or an
  array, and test one call of any watcher's probe before backgrounding it.
- **Pytest in the container:** the container has no pytest of its own (never pip-install into it; use a
  scratch `pylib`). Always pass `--basetemp` and set `TMPDIR` under `/data/.../mh_scratch`, and remove them
  afterwards: the container's `/tmp` is on appct's small `/var`, which full suites have filled (Errno 28).

## Working with Matthew

- Merge, push tags and deploy only on an explicit request for that specific action.
- Report outcomes as measured: a "done" job with exit 0 can still have lost frames, so check its log.
- Background reprocessing (P2/P3) can wait; nightly data (P0/P1) cannot.
