PYTHON=python3
INC=-Ilib -I${CONDA_PREFIX}/include
EXTRA=-O3
VERSION=$(shell git show -s --pretty="%h %ad %d")
CFLAGS=-fPIC -std=c++11 -Wall -g -fopenmp ${INC}
SO=_sib$(shell ${PYTHON}-config --extension-suffix)
LINK=-lgomp -lm -DVERSION="\"${VERSION}\""
PYINC=$(shell ${PYTHON} -m pybind11 --includes)
DEP=$(wildcard *.h)
DOCTEST=$(wildcard test/*.doctest)
CXX=g++

all: sib ${SO}

params.o: params.cpp ${DEP}
	${CXX} ${CFLAGS} ${EXTRA} -c params.cpp -o $@
bp.o: bp.cpp ${DEP}
	${CXX} ${CFLAGS} ${EXTRA} -c bp.cpp -o $@
sib: bp.o params.o sib.cpp ${DEP}
	${CXX} ${CFLAGS} ${EXTRA} params.o bp.o sib.cpp ${LINK} -o $@
drop.o: drop.cpp ${DEP}
	${CXX} ${CFLAGS} ${EXTRA} -c drop.cpp -o $@
${SO}: bp.o params.o drop.o pysib.cpp ${DEP}
	${CXX}  -shared ${CFLAGS} ${PYINC} ${LINK} ${EXTRA} params.o bp.o drop.o pysib.cpp -o $@

test: all doctest
	${PYTHON} test/run_tests.py
doctest:
	for x in ${DOCTEST}; do ${PYTHON} -c "import sys, doctest; (f,t) = doctest.testfile(\"$$x\"); print(f'DOCTEST $$x: PASSED {t-f}/{t}'); sys.exit(int(f > 0))"; done

clean:
	rm -f sib ${SO} *.o

.PHONY: test
