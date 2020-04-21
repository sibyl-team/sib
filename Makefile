CFLAGS=-Wall -O3 -g

sib: sib.cpp factorgraph.cpp factorgraph.h
	g++ ${CFLAGS} factorgraph.cpp sib.cpp -lm -o $@
clean:
	rm -f sib
