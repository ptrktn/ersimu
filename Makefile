OCTAVE = octave

# additional compiler flags
CFLAGS += -O3

# silence certain warnings in lsoda.c
LSODACFLAGS=-Wno-maybe-uninitialized -Wno-unused-but-set-variable \
	-Wno-unused-variable

.PHONY: test
test: test-lv

.PHONY: test-lv
test-lv:
	rm -f lotkavolterra.pdf
	./ersimu.py --verbose --scipy --name lotkavolterra --run examples/lotka-volterra.txt
	test -f lotkavolterra.pdf


.PHONY: test-2
test-2:
	rm -f ersimu.mat ersimu.m ersimu_scipy.py ersimu.pdf
	./ersimu.py -v examples/szalai-koros-chd.txt
	./ersimu_scipy.py
	test -f ersimu.pdf

.PHONY: test-o
test-o:
	rm -f ersimu.mat ersimu.m ersimu.h
	./ersimu.py examples/oregonator.txt
	$(OCTAVE) ersimu.m
	./xplot.sh -N test-o ersimu.mat
	rm -f ersimu.mat
	$(MAKE) lsoda
	./lsoda > ersimu.mat
	./xplot.sh ersimu.mat

lsoda: lsoda.c ersimu.h
	$(CC) $(CFLAGS) $(LSODACFLAGS) -o lsoda lsoda.c -lm

.PHONY: test-lsoda
test-lsoda:
	rm -f ersimu.h ersimu.mat
	./ersimu.py examples/brusselator.txt
	$(MAKE) lsoda
	./lsoda > ersimu.mat
	./xplot.sh -N test-lsoda ersimu.mat

.PHONY: lsoda-sim
lsoda-sim: lsoda
	rm -f ersimu.mat
	./lsoda > ersimu.mat

.PHONY: clean
clean:
	rm -f core *~ *.BAK octave-workspace *.jpg *.mat *.m *.h lsoda

.PHONY: test-parser
test-parser:
	@for i in $(shell find examples -type f -name '*.txt' | sort) ; do \
         echo "Input file: $$i" ; \
         ./ersimu.py $$i || exit 1 ; \
    done

.PHONY: dep
dep:
	pip3 install -r requirements.txt --user
