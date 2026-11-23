# Development

For building and running the adjoint, see the [README](../README.md).

## Debugging the build

The Dockerfile is a multi-stage build. The `builder` stage installs
`openmethane-cmaq-dev` into a clean `debian:trixie-slim` and compiles the
adjoint binaries, but it isn't part of any exported target, so a compile
failure leaves nothing to look at. To debug it, use:

```shell
# build only the "builder" stage and tag as "cmaq-adjoint-builder"
docker build . --progress=plain --target builder -t cmaq-adjoint-builder

# run an interactive bash terminal in the builder
docker run -it --rm cmaq-adjoint-builder bash
```

The sources are in `/opt/cmaq-adjoint/cmaq`, the build directories in
`/opt/cmaq-adjoint/scripts`, and the binaries in `/opt/cmaq/bin`. To build
against an unmerged `openmethane-cmaq-dev`, add
`--build-context cmaq-dev=<directory holding the .deb>`. The steps are
described in [building.md](building.md).

## Testing

A simple test script can be invoked with:

```shell
make test
```

This builds the `test` stage, a clean `debian:trixie-slim` with only the
package installed, which has already checked with `ldd` that every ELF file in
the package finds its libraries. It then mounts the `tests` folder under
`/opt/tests` and runs `/opt/tests/test-run.sh`, as CI does. The test script
downloads a bundle of test data and a file describing necessary environment
variables, and then runs the forward and backward adjoint on this input.

**Note:** this does **not** perform a numerical test, it simply ensures that
the compiled binary will successfully run to completion on known inputs.

## Code review

Review isn't only for catching bugs. It keeps the team aware of how the
project is changing and builds agreement on decisions. Not every pull request
benefits from that. The author of a pull request decides whether it needs a
human review, and says so on the first line of the description, where it
shows up in email notifications:

```markdown
**Review:** needed — changes science outputs: the adjoint forcing is now read at the end of the output step
```

```markdown
**Review:** not needed — test harness only, no model source changes; merging Thursday unless someone objects
```

Always give a reason. The reason lets everyone else check the decision at a
glance, and shows over time whether these categories need adjusting.

The lists below are a guide, not rules, and some judgement is still needed.
When a change doesn't clearly fit, ask: would the rest of the team want to
know about this, or have a say in it? If so, it needs review.

### Review needed

- **Changes to science outputs.** Any change to what `ADJOINT_FWD` or
  `ADJOINT_BWD` produces needs review, however small the fix that causes it.
  That covers the forward concentrations, the adjoint sensitivities and
  gradients, and anything that changes the numerics on the way to them, such
  as checkpointing, the timing of the forcing or emissions, or the headers of
  output files that downstream code reads.
- **Changes to the build that can change results**, such as a new base image,
  compiler or compiler flags. A smoke test passing is not evidence that the
  numbers are unchanged.
- **Refactors that change structure or behaviour across more than one file.**
- **New concepts or features**, including new outputs, environment variables
  and methodology decisions.
- **Documentation changes.** If something is worth documenting, it's worth
  the rest of the team knowing about it and agreeing with it. This doesn't
  include fixes to typos, grammar or formatting, which are mechanical changes.
- **Changes that make a production run fail, or stop failing, where it didn't
  before.** Make sure whoever looks after production runs knows about these.

### Review not needed

- **Infrastructure changes** made by the person who looks after the
  infrastructure.
- **Mechanical changes**, such as small corrections and bug fixes whose
  correctness is established by tests, tool-generated changes such as
  running a formatter, and fixes to typos, grammar or formatting in
  documentation. Being mechanical is what matters here, not size: a
  formatter run that touches every file is still mechanical.
- **Changes to the test harness** that don't touch the model source or build.
- **Small changes to logging, warnings and project settings.**

### Science code isn't the same as science outputs

A change can touch model source without changing any outputs. For example,
it can add a check that stops a misconfigured run, while a correctly
configured run produces the same files as before. These changes don't need
review, but the reason has to give the evidence that the outputs don't
change, not just say so. If you can't write that sentence, choose "needed".

### Agent review before a pull request

If you use an agent such as Claude Code to prepare a change, `CLAUDE.md`
asks it to review anything that isn't mechanical with a fresh reviewer
(`.claude/agents/reviewer.md`) before it opens the pull request. The agent
fixes clear problems itself and brings anything else to you. For riskier
changes it offers a second, independent review. Each review costs tokens, so
you can tell the agent to skip it. This is part of preparing the change, not
a substitute for human review: it doesn't change whether a pull request needs
review.

### Not needed doesn't mean merged immediately

For a "not needed" pull request, say when you'll merge it. This gives
everyone a chance to ask for a review before then. Anyone can switch a pull
request to "needed", and then it waits for a review as usual.

### Examples

| Pull request | First line |
|---|---|
| [#20](https://github.com/openmethane/CMAQ-Adjoint/pull/20) | **Review:** not needed — adds a backward smoke test; only `tests/` changes, no model source or build changes |
| [#23](https://github.com/openmethane/CMAQ-Adjoint/pull/23) | **Review:** not needed — restores a check that the backward run's emissions layer count matches the emissions file; the production configuration passes the smoke test unchanged |
| [#19](https://github.com/openmethane/CMAQ-Adjoint/pull/19) | **Review:** needed — changes science outputs: the adjoint forcing is now read at the end of the output step, one step earlier than before |
| [#28](https://github.com/openmethane/CMAQ-Adjoint/pull/28) | **Review:** needed — build change that can change results: pins the base image, which moves the compiler from gfortran 12 to 14.2 |
| [#33](https://github.com/openmethane/CMAQ-Adjoint/pull/33) | **Review:** needed — new feature: the backward run can write the sensitivity to the lateral boundary condition, plus docs |
| [#29](https://github.com/openmethane/CMAQ-Adjoint/pull/29) | **Review:** needed — documentation: findings on how linear the CH4 model is |

#23 is "not needed" because a correctly configured run is unaffected.
Because a misconfigured run now stops where it used to run to completion, the
person who looks after production runs should still hear about it.
