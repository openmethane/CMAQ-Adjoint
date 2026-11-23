# Debian 13 (trixie) on amd64. See docs/building.md for the native build.
TEST_IMAGE=cmaq-adjoint-test

.phony: build
build: Dockerfile ## build openmethane-cmaq-adjoint, exported to dist/
	docker build --platform=linux/amd64 --progress=plain --target deb --output type=local,dest=dist .

.phony: test
test: Dockerfile ## install the package in a clean image and run tests/test-run.sh
	docker build --platform=linux/amd64 --progress=plain --target test -t ${TEST_IMAGE} .
	docker run --rm -v ./tests:/opt/tests ${TEST_IMAGE} /opt/tests/test-run.sh
