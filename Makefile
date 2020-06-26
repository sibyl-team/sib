INC=-Ilib $(patsubst %, -I%/include, ${CONDA_PREFIX})
VERSION=$(shell git show -s --pretty="%h %ad %d")
CFLAGS=-fPIC -std=c++11 -Wall -O3 -g -fopenmp ${INC}
SO=build/_sib$(shell python3-config --extension-suffix)
LINK=-lgomp -lm -DVERSION="\"${VERSION}\""
PYINC=$(shell python3 -m pybind11 --includes)
BUILD=build
OBJS=$(patsubst %, ${BUILD}/%, params.o bp.o drop.o)
CXX=g++

all: build sib ${SO}
build:
	mkdir -p build
${OBJS}: build/%.o : src/%.cpp src/params.h
	${CXX} ${CFLAGS} -c $< -o $@
sib: ${OBJS} build src/sib.cpp
	${CXX} ${CFLAGS} ${OBJS} src/sib.cpp ${LINK} -o $@
${SO}: ${OBJS} build src/pysib.cpp
	${CXX}  -shared ${CFLAGS} ${PYINC} ${LINK} ${OBJS} src/pysib.cpp -o $@

test: all
	python3 test/run_tests.py


clean:
	rm -f sib ${OBJS} ${SO}

.PHONY: test
