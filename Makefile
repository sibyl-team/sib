CFLAGS=-Wall -O3 -g

sib: sib.cpp
	g++ ${CFLAGS} sib.cpp -lm -o $@
clean:
	rm -f sib
