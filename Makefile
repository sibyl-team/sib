CFLAGS=-Wall -O3 -g

sib: sib.cpp bp.cpp bp.h params.o cavity.h
	g++ ${CFLAGS} bp.cpp sib.cpp params.o -lm -o $@
clean:
	rm -f sib
