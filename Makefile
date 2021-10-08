OCTAVE = octave
SHELL = bash

# silence certain warnings in lsoda.c
LSODACFLAGS = -Wno-maybe-uninitialized -Wno-unused-but-set-variable \
	-Wno-unused-variable

# additional compiler flags
CFLAGS += -O3 $(LSODACFLAGS)

.PHONY: test
test:
	$(MAKE) PYTHONDONTWRITEBYTECODE=1 test-all

.PHONY: test-all
test-all: test-parse test-lsodac test-scipy test-octave

.PHONY: test-parse
test-parse:
	for i in `find examples -type f | sort` ; do \
      rm -f ersimu.h simulation.m simulation.py simulation.tex ; \
      for opt in --lsodac --octave --scipy --latex ; do \
          echo Parse $$i option $$opt ; \
          ./ersimu.py $$opt $$i || exit 1 ; \
          name=`basename $$i` ; \
          name=$${name%.*} ; \
          name=$${name//-/} ; \
          echo Parse $$i option $$opt --name $$name ; \
          ./ersimu.py $$opt $$i --name $$name || exit 1 ; \
      done ; \
      pdflatex $$name.tex > /dev/null 2>&1 || exit 1 ; \
      pdflatex $$name.tex > /dev/null 2>&1 || exit 1 ; \
      test -f ersimu.h || exit 1 ; \
      test -f simulation.m || exit 1 ; \
      test -f simulation.py || exit 1 ; \
      test -f simulation.tex || exit 1 ; \
      pdflatex simulation.tex > /dev/null 2>&1 || exit 1 ; \
    done

.PHONY: test-scipy
test-scipy:
	rm -f test.pdf
	./ersimu.py --verbose --scipy --name test --run examples/oregonator.txt
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
	test -f simulation.dat
	./xplot.sh -N test simulation.dat

%.o: %.c
	$(CC) $(CFLAGS) -c $<

simulation: simulation.o
	$(CC) $? -o $@ -lm

.PHONY: test-lsodac
test-lsodac:
	./ersimu.py examples/brusselator.txt
	$(MAKE) simulation
	./simulation test.dat
	./xplot.sh -N test test.dat

.PHONY: clean
clean:
	rm -f core *~ *.o *.BAK *.dat *.m ersimu.c *.pdf simulation

.PHONY: dep
dep:
	pip3 install -r requirements.txt --user

.PHONY: install-platform
install-platform:
	sudo apt install build-essential gnuplot python3 python3-pip \
         texlive-latex-base texlive-latex-extra \
         octave octave-doc
