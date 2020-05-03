CFLAGS=-fPIC -std=c++11 -Wall -O3 -g -fopenmp
SO=_sib$(shell python3-config --extension-suffix)
LINK=-lgomp -lm
PYINC=$(shell python3 -m pybind11 --includes)
CXX=g++

all: sib ${SO}

sib: sib.cpp bp.h bp_impl.h params.h params_impl.h
	${CXX} ${CFLAGS} sib.cpp ${LINK} -o $@
${SO}: pysib.cpp bp.h bp_impl.h params.h params_impl.h
	${CXX}  -shared ${CFLAGS} ${PYINC} ${LINK} pysib.cpp -o $@

test: all
	python3 test/run_tests.py


clean:
	rm -f sib ${SO} *.o

.PHONY: test
