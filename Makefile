CFLAGS=-Wall -O3 -g

sib: sib.cpp bp.cpp bp.h
	g++ ${CFLAGS} bp.cpp sib.cpp -lm -o $@
clean:
	rm -f sib
