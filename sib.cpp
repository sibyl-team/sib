
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
#include "factorgraph.h"

using namespace std;

int main(int argc, char ** argv) {
	FactorGraph factor(Params(argc, argv));
	return 0;
}


