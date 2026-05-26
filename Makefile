IMAGE=cmaq-adjoint

.phony: build
build: Dockerfile ## build an x86_64 version of the docker container
	docker build --platform=linux/amd64 --progress=plain . -t ${IMAGE}

.phony: build-aarch64
build-aarch64: Dockerfile  ## build an arm version of the docker container
	docker build --platform=linux/arm64 . -t ${IMAGE}

.phony: run
run: build  ## run the docker container
	docker run -it --rm -v .:/opt/project ${IMAGE}

.phony: test
test: build
	docker run -it --rm -v ./tests:/opt/tests ${IMAGE} /opt/tests/test-run.sh
