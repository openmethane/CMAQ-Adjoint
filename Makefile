.phony: build
build: Dockerfile ## build an x86_64 version of the docker container
	docker build --platform=linux/amd64 --progress=plain . -t cmaq

.phony: build-aarch64
build-aarch64: Dockerfile  ## build an arm version of the docker container
	docker build --platform=linux/arm64 . -t cmaq

.phony: build-conda
build-conda: Dockerfile  ## build just the "conda" step of the docker container
	docker build --platform=linux/arm64 . --target conda -t cmaq-conda

.phony: run
run: build  ## run the docker container
	docker run -it --rm -v ${PWD}:/opt/project cmaq

.phony: run-conda
run-conda: build-conda  ## run a container with only the conda dependencies
	docker run -it --rm -v ${PWD}:/opt/project cmaq-conda

.phone: test
test: build
	docker run -it --rm -v ${PWD}/tests:/opt/tests cmaq /opt/tests/test-run.sh
