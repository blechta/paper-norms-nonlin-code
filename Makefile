DOCKER_IMG=quay.io/blechta/fenics-dev:paper0
DOCKER_INIT=docker create -v /home/fenics/.instant --name instant-cache \
	       $(DOCKER_IMG) /bin/true 2>/dev/null || true
DOCKER_RUN=docker run --volumes-from instant-cache --rm \
		   -v $(shell pwd):/home/fenics/work -w /home/fenics/work \
		   $(DOCKER_IMG) "sudo /bin/bash -l -c '$(CMD)'"

INIT=$(DOCKER_INIT)
RUNNER=$(DOCKER_RUN)

LOGS= \
    CarstensenKlose_4.0_05.log \
    CarstensenKlose_4.0_10.log \
    CarstensenKlose_4.0_15.log \
    CarstensenKlose_4.0_20.log \
    ChaillouSuri_10.0_05.log \
    ChaillouSuri_10.0_10.log \
    ChaillouSuri_10.0_15.log \
    ChaillouSuri_10.0_20.log \
    ChaillouSuri_1.5_05.log \
    ChaillouSuri_1.5_10.log \
    ChaillouSuri_1.5_15.log \
    ChaillouSuri_1.5_20.log

.PHONY: all

all: tabular.tex

tabular.tex: parse_results.py $(LOGS)
	python $< > $@

%.log: main.py
	$(eval CMD=python $< $(subst _, ,$*) 2>&1 > $@ || true)
	$(INIT)
	$(RUNNER)

clean-all: clean-local clean-docker

clean-local:
	rm -f $(LOGS) tabular.tex

clean-cache:
	docker rm instant-cache || true

clean-docker: clean-cache
	docker rmi $(DOCKER_IMG)
