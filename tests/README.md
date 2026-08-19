# Testing

These tests check whether the adjoint binaries compile and run to completion
on known inputs. They are smoke tests: they do not verify that the numbers
coming out are correct.

`tests/test-run.sh` runs three cases:

- **001** — `ADJOINT_FWD`, single process. Checks normal completion.
- **002** — `ADJOINT_FWD` under `mpirun -np 2`. Checks normal completion.
- **003** — `ADJOINT_FWD` then `ADJOINT_BWD`, single process. Checks normal
  completion, no non-finite values in the log, and that the sensitivity output
  is not empty.

The backward run consumes the forward run's checkpoints, `XFIRST` state and
concentration output, so case 003 runs both binaries in the same directory
rather than starting the backward run from clean.

## Environment

Both env files ship in the test-data bundle:

- `fwd.env` — the domain, time period and input paths. Used by all three
  cases. (Previously `test.env`.)
- `bwd.env` — everything specific to the backward run. Case 003 sources it
  *after* `fwd.env`, so it inherits the shared settings and overrides only
  what differs.

`bwd.env` mirrors py4dvar's `run_bwd_single()` in
`openmethane/src/openmethane/fourdvar/util/cmaq_handle.py`, which is what
drives this binary in production. Keep the two in step.

Two settings are deliberately *not* inherited from `fwd.env`, both because
`FLCHECK` rejects any input file whose start time is absent:

- `INIT_*_1` stay pointed at the initial conditions rather than the forward
  run's final `CGRID`. py4dvar uses the previous day's `CGRID`, which in its
  multi-day chain is stamped at the current day's start; a single-day test
  stamps it at the run end.
- `S_CGRID` is redirected to a backward-specific path, since the forward run
  has just written that file into the shared directory.

Because the env files live in the bundle, changing the backward configuration
means rebuilding and re-uploading the bundle rather than editing the repo.

## The adjoint forcing file

Case 003 needs an `ADJ_FORCE` input. `FORCE.20221207.nc` in the test-data
bundle is all zeros except a single cell in its **final** time step.

That shape is deliberate. The backward driver walks output steps from the run
end towards the run start, so the final time step is the first one it reads. A
forcing spread across every time step would leave the test unable to distinguish
a correct run from one that reads the time step one output step out of step —
both would produce a plausible non-zero field. With one populated time step at
the run end, an off-by-one read misses it completely and the sensitivities
come back identically zero.

Measured on this domain: the correct run writes ~245 kB of non-zero bytes into
`ADJ_LGRID`, a run that never picks the forcing up writes ~6 kB (the netCDF
header alone). Note that the latter still reports *normal completion* — the
emptiness check, not the completion check, is what catches it.

To regenerate the file:

```shell
# needs netCDF4, which the runtime image deliberately does not carry;
# py4dvar's virtualenv has it
/path/to/openmethane/.venv/bin/python tests/make-force-file.py \
    --template /path/to/openmethane/tests/test-data/templates/force_template.nc \
    --output tests/test-data/FORCE.20221207.nc
```

py4dvar's `force_template.nc` is used as the base because it is the same
domain (10x10x32, `GDNAM` `au-test_v1`), the same date, and the same 25 hourly
records as the rest of this bundle — i.e. exactly what the model is fed in
production.

## Test data

As the test data is too large to be committed to the repo, it is fetched from
the Open Methane Public Data Store. A cached copy in `test-data/` is reused if
present, so delete that directory to pick up a new bundle.

The bundle contains:

```
fwd.env                 environment shared by all cases
bwd.env                 backward-run environment, sourced after fwd.env
FORCE.20221207.nc       adjoint forcing (see above)
mcip/                   sample MCIP output for the au-test_v1 domain
emissions/              emissions for the matching domain and date
cmaq/                   initial and boundary conditions
```

Test data is bundled by compressing the contents of `test-data/`. Write the
archive outside the repo so it is not picked up as an untracked file, and
date-stamp it so that a change to the bundle does not silently alter what
older commits fetch:

```shell
# from tests/
tar -C test-data -czf ../../cmaq-adjoint-test-data.$(date +%Y%m%d).tar.gz .
```

Then upload it to the Open Methane Public Data Store and point
`TEST_DATA_URL` in `test-run.sh` at the new file. Because the URL is
versioned, updating the bundle is a code change: the two land together.

`test-run.sh` checks for `fwd.env`, `bwd.env` and `FORCE.20221207.nc` up
front, so a cached bundle predating the backward test fails with a clear
message rather than part way through a model run.
