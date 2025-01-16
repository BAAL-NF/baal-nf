LOCAL_IMAGE := baal-nf-env
REMOTE_IMAGE := oalmelid/${LOCAL_IMAGE}
TAG := latest
.DEFAULT_GOAL := image

.PHONY: image publish

image: .image_built

.image_built: Dockerfile env/* 
	DOCKER_BUILDKIT=1 docker build -t ${LOCAL_IMAGE}:${TAG} . &&\
	touch .image_built

publish: image
	docker tag $(LOCAL_IMAGE):$(TAG) $(REMOTE_IMAGE):$(TAG)
	docker push $(REMOTE_IMAGE):$(TAG)

