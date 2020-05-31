INC=-I${CONDA_PREFIX}/include
CFLAGS=-fPIC -std=c++11 -Wall -O3 -g -fopenmp ${INC}
SO=_sib$(shell python3-config --extension-suffix)
LINK=-lgomp -lm
PYINC=$(shell python3 -m pybind11 --includes)
CXX=g++

all: sib ${SO}

params.o: params.cpp params.h
	${CXX} ${CFLAGS} -c params.cpp -o $@
bp.o: bp.cpp bp.h cavity.h
	${CXX} ${CFLAGS} -c bp.cpp -o $@
sib: bp.o params.o sib.cpp 
	${CXX} ${CFLAGS} params.o bp.o sib.cpp ${LINK} -o $@
drop.o: drop.cpp
	${CXX} ${CFLAGS} -c drop.cpp ${LINK} -o $@
${SO}: bp.o params.o drop.o pysib.cpp
	${CXX}  -shared ${CFLAGS} ${PYINC} ${LINK} params.o bp.o drop.o pysib.cpp -o $@

test: all
	python3 test/run_tests.py


clean:
	rm -f sib ${SO} *.o

.PHONY: test
