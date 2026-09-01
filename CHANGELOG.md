# Changelog

Versions follow [Semantic Versioning](https://semver.org/) (`<major>.<minor>.<patch>`).

Backward incompatible (breaking) changes will only be introduced in major versions
with advance notice in the **Deprecations** section of releases.


<!--
You should *NOT* be adding new changelog entries to this file, this
file is managed by towncrier. See changelog/README.md.

You *may* edit previous changelogs to fix problems like typo corrections or such.
To add a new changelog entry, please see
https://pip.pypa.io/en/latest/development/contributing/#news-entries,
noting that we use the `changelog` directory instead of news, md instead
of rst and use slightly different categories.
-->

<!-- towncrier release notes start -->

## CMAQ Adjoint v2.0.1 (2026-09-01)

### 🎉 Improvements

- Pin the Docker base image to `ghcr.io/openmethane/cmaq:1.0.2` (Debian trixie) instead of the floating `stable` tag, and add `-fbacktrace` to the forward build's link flags to match the backward build. Note that the trixie base compiles with gfortran 14 rather than 12, so numerical output may differ slightly from previous builds. ([#28](https://github.com/openmethane/CMAQ-adjoint/issues/28))


## CMAQ Adjoint v2.0.0 (2026-08-19)

### 🎉 Improvements

- Include zenodo versions of science modules ([#14](https://github.com/openmethane/CMAQ-adjoint/issues/14))
- Added a backward (`ADJOINT_BWD`) smoke test, which runs the backward adjoint against a forward run's checkpoints and a small adjoint forcing file, and checks that it completes and produces non-empty sensitivities. The test data bundle gained a forcing file and a backward environment file, and `test.env` was renamed to `fwd.env`. ([#20](https://github.com/openmethane/CMAQ-adjoint/issues/20))

### 🐛 Bug Fixes

- Corrected the `ADJ_LGRID_EM` emissions sensitivity file header, which described a different layer count and species list than the array actually written to it. The header is now built from the same rule the backward vertical diffusion uses to fill that array, and the backward run stops with an error if the two ever disagree. The `units` attribute is now reported as `CF/(mol/s)`, matching the emissions file units the sensitivity is taken with respect to. ([#18](https://github.com/openmethane/CMAQ-adjoint/issues/18))
- Corrected the time at which the backward run reads each adjoint forcing record. The record is now read at the end of the output step it is applied to, matching the times at which the forward run writes concentrations, rather than one output step earlier. This changes computed sensitivities. ([#19](https://github.com/openmethane/CMAQ-adjoint/issues/19))
- Stopped the `CREATE_CHK` environment variable, which suppresses checkpoint creation for finite-difference testing in the forward run, from also suppressing the backward run's opening of those checkpoint files for reading. ([#22](https://github.com/openmethane/CMAQ-adjoint/issues/22))
- The backward run now checks that the number of emissions layers it was configured with, via `CTM_EMLAYS`, matches the emissions file, and stops with an error if it does not. The forward run takes this value from the emissions file directly, so a mismatched setting previously gave the adjoint a different number of emissions layers than the forward run used, without reporting anything. ([#23](https://github.com/openmethane/CMAQ-adjoint/issues/23))
- fix emission interpolation and sensitivity accumulation to match py4dvar ([#25](https://github.com/openmethane/CMAQ-adjoint/issues/25))
- combine continuous and discrete adjoint to produce correct gradient ([#27](https://github.com/openmethane/CMAQ-adjoint/issues/27))


## CMAQ Adjoint v1.0.1 (2026-05-26)

No significant changes.


## CMAQ Adjoint v1.0.0 (2026-05-26)

### 🆕 Features

- Initial fork of CMAQ Adjoint from the [development repo](https://adjoint.colorado.edu:8080) at University of
  Colorado Boulder

### 🎉 Improvements

- Add simple test script ([#2](https://github.com/openmethane/CMAQ-adjoint/issues/2))
- Build and push docker image in GitHub Actions ([#2](https://github.com/openmethane/CMAQ-adjoint/issues/2))
- Build CMAQ libraries and dependencies in docker ([#2](https://github.com/openmethane/CMAQ-adjoint/issues/2))
- Speed up tests by using the smaller au-test domain ([#5](https://github.com/openmethane/CMAQ-adjoint/issues/5))
- Use [`openmethane/cmaq`](https://github.com/openmethane/CMAQ) as the base docker image, which provides:
  - build dependencies such as mpich and libnetcdf
  - libraries such as ioapi, pario and stenex
  - build configuration for CMAQ bldit scripts
  - binaries for MCIP, BCON_CH4only and ICON CH4only

  ([#9](https://github.com/openmethane/CMAQ-adjoint/issues/9))
