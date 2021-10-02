OCTAVE = octave

# additional compiler flags
CFLAGS += -O3

# silence certain warnings in lsoda.c
LSODACFLAGS=-Wno-maybe-uninitialized -Wno-unused-but-set-variable \
	-Wno-unused-variable

.PHONY: test
test:
	$(MAKE) PYTHONDONTWRITEBYTECODE=1 test-all

.PHONY: test-all
test-all: test-parse test-lsodac test-scipy test-octave

.PHOHY: test-parse
test-parse:
	for i in `find examples -type f | sort` ; do \
      for opt in --lsodac --octave --scipy --latex ; do \
          echo Parse $$i option $$opt ; \
          ./ersimu.py $$opt $$i || exit 1 ; \
      done ; \
    done

.PHONY: test-scipy
test-scipy:
	rm -f test.pdf
	./ersimu.py --verbose --scipy --name test --run examples/lotka-volterra.txt
	test -f test.pdf

.PHONY: test-2
test-2:
	rm -f ersimu.mat ersimu.m ersimu_scipy.py ersimu.pdf
	./ersimu.py -v examples/szalai-koros-chd.txt
	./ersimu_scipy.py
	test -f ersimu.pdf

.PHONY: test-octave
test-octave:
	rm -f simulation.dat simulation.m
	./ersimu.py --octave examples/oregonator.txt
	test -f simulation.m
	$(OCTAVE) simulation.m
	./xplot.sh -N test simulation.dat

lsoda: lsoda.c ersimu.h
	$(CC) $(CFLAGS) $(LSODACFLAGS) -o lsoda lsoda.c -lm

.PHONY: test-lsodac
test-lsodac:
	rm -f ersimu.h ersimu.mat
	./ersimu.py examples/brusselator.txt
	$(MAKE) lsoda
	./lsoda test.dat
	./xplot.sh -N test test.dat

.PHONY: clean
clean:
	rm -f core *~ *.BAK octave-workspace *.jpg *.mat *.m *.h *.pdf lsoda

.PHONY: dep
dep:
	pip3 install -r requirements.txt --user
