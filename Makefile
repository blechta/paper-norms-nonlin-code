DOCKER_IMG=quay.io/blechta/dolfin-tape@sha256:ee02b293f1e265011e7abe07cf53b1808fc1ae0cef603f26915db4fa4c9b16bc
DOCKER_CACHE=instant-cache
DOCKER_RUN=docker run --volumes-from $(DOCKER_CACHE) --rm \
		   -v $(shell pwd):/home/fenics/work -w /home/fenics/work \
		   $(DOCKER_IMG) "sudo /bin/bash -l -c '$(CMD)'"

INIT=init-docker
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

.PHONY: all init-docker clean-all clean-local clean-cache clean-docker

.PRECIOUS: %.log

all: tabular.tex

tabular.tex: parse_results.py $(LOGS)
	python $< > $@

%.log: main.py $(INIT)
	$(eval CMD=python $< $(subst _, ,$*) >$@ 2>&1 || true)
	$(RUNNER)

init-docker:
	docker history -q $(DOCKER_IMG) >/dev/null 2>&1 || docker pull $(DOCKER_IMG)
	docker inspect $(DOCKER_CACHE) >/dev/null 2>&1 || \
	    docker create -v /home/fenics/.instant --name $(DOCKER_CACHE) $(DOCKER_IMG) /bin/true

clean-all: clean-local clean-docker

clean-local:
	rm -f $(LOGS) *.pdf tabular.tex

clean-cache:
	docker rm $(DOCKER_CACHE) >/dev/null 2>&1 || true

clean-docker: clean-cache
	docker rmi $(DOCKER_IMG) >/dev/null 2>&1 || true
