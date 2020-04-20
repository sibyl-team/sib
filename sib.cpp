
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <iostream>
#include "siblib.h"

using namespace std;

int main(int argc, char ** argv) {

	Params par;
	par.read_input(argc, argv);

	Factor factor(&par);

	return 0;
}


