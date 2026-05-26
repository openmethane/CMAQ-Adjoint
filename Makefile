.phony: build
build: Dockerfile ## build an x86_64 version of the docker container
	docker build --platform=linux/amd64 --progress=plain . -t cmaq-adj

.phony: build-aarch64
build-aarch64: Dockerfile  ## build an arm version of the docker container
	docker build --platform=linux/arm64 . -t cmaq-adj

.phony: run
run: build  ## run the docker container
	docker run -it --rm -v ${PWD}:/opt/project cmaq-adj

.phony: test
test: build
	docker run -it --rm -v ${PWD}/tests:/opt/tests cmaq-adj /opt/tests/test-run.sh
