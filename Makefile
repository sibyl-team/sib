CFLAGS=-Wall -fPIC -O3 -g
SO=sib$(shell python3-config --extension-suffix)

all: sib ${SO}

bp.o: bp.cpp bp.h params.o cavity.h
	c++ ${CFLAGS} -c bp.cpp -o $@
sib: bp.o params.o
	c++ bp.o params.o sib.cpp -lm -o sib
${SO}: bp.o
	c++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` pysib.cpp bp.cpp params.cpp -o ${SO}

clean:
	rm -f sib ${SO} *.o
