.phony: build
build: Dockerfile ## build an x86_64 version of the docker container
	docker build --platform=linux/amd64 . -t cmaq

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


tests/test-data/mcip:  # The required MCIP output is too large to store in the repository so fetch it from a S3 bucket
	wget -c -nv https://prior.openmethane.org/cmaq/tests/test-data/mcip.tar.gz -O tests/test-data/mcip.tar.gz
	# Extracts to tests/test-data/mcip
	tar -xf tests/test-data/mcip.tar.gz


upload-mcip:  # Update the MCIP data in the S3 bucket
	tar -czvf tests/test-data/mcip.tar.gz tests/test-data/mcip
	aws s3 cp tests/test-data/mcip.tar.gz s3://openmethane-prior/cmaq/tests/test-data/mcip.tar.gz --endpoint-url https://8f8a25e8db38811ac9f26a347158f296.r2.cloudflarestorage.com --profile cf-om-prior-r2
