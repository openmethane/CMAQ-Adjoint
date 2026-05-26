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
