LOCAL_IMAGE := baal-chip-env
REMOTE_IMAGE := oalmelid/baal-chip-env
TAG := latest
RANDOM := $(shell bash -c 'echo $$RANDOM')
.DEFAULT_GOAL := image
.PHONY: image publish

image: .docker.log

.docker.log: Dockerfile env/* scripts/*
	DOCKER_BUILDKIT=1 docker build --build-arg BUILD_FORCE=${RANDOM} --progress=plain --secret id=gitlab,src=credentials.sh -t baal-chip-env:$(TAG) . 2>&1 | tee .docker.log

publish: image
	docker tag $(LOCAL_IMAGE):$(TAG) $(REMOTE_IMAGE):$(TAG)
	docker push $(REMOTE_IMAGE):$(TAG)
 
