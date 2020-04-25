CFLAGS=-fPIC -std=c++11 -Wall -O3 -g -fopenmp
SO=_sib$(shell python3-config --extension-suffix)
LINK=-lgomp -lm
PYINC=$(shell python3 -m pybind11 --includes)
CXX=g++

all: sib ${SO}

bp.o: bp.cpp bp.h cavity.h
	${CXX} ${CFLAGS} -c bp.cpp -o $@
sib: bp.o sib.cpp
	${CXX} ${CFLAGS} bp.o sib.cpp ${LINK} -o $@
${SO}: bp.o pysib.cpp
	${CXX}  -shared ${CFLAGS} ${PYINC} ${LINK} pysib.cpp bp.o -o $@

clean:
	rm -f sib ${SO} *.o
