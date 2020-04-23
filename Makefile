CFLAGS=-Wall -fPIC -O3 -g

all: sib sib.so

bp.o: bp.cpp bp.h params.o cavity.h
	c++ ${CFLAGS} -c bp.cpp -o $@
sib: bp.o params.o
	c++ bp.o params.o sib.cpp -lm -o sib
sib.so: bp.o
	c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` pysib.cpp bp.cpp params.cpp -o sib`python3-config --extension-suffix`

clean:
	rm -f sib
