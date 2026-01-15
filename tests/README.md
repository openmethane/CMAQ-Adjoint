# Testing

These tests only check whether the forward adjoint has compiled and is able
to run without errors on known inputs.

## Test data

The test data provided is sample output from MCIP, plus an initial condition
and boundary condition for the matching domain.

As the test data is too large to be committed to the repo, it is fetched from
the Open Methane Public Data Store.

Test data is bundled using the command:

```shell
# compress the contents of test-data/ into cmaq-adj-test-data.tar.gz
tar -C test-data -czf cmaq-adj-test-data.tar.gz .
```
