CFLAGS=-std=c++11 -Wall -O3 -g -fopenmp -lgomp
SO=_sib$(shell python3-config --extension-suffix)
PYINC=$(shell python3 -m pybind11 --includes)

all: sib ${SO}

bp.o: bp.cpp bp.h params.o cavity.h
	c++ ${CFLAGS} -c bp.cpp -o $@
sib: bp.o params.o
	c++ ${CFLAGS} bp.o params.o sib.cpp -lm -o sib
${SO}: bp.cpp pysib.cpp params.cpp
	c++ -shared -fPIC ${CFLAGS} ${PYINC} pysib.cpp bp.cpp params.cpp -o ${SO}

clean:
	rm -f sib ${SO} *.o
