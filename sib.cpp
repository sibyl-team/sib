
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
#include "bp.h"

using namespace std;

int main(int argc, char ** argv) {
	FactorGraph factor(Params(argc, argv));
	factor.init_msg();
	factor.showmsg();
	int N = int(factor.nodes.size());
	real_t err = 1.0;
	int it = 0;
	while(err > factor.params.tol) {
		++it;
		err = 0.0;
		for(int i = 0; i < N; ++i) {
			real_t err_i = factor.update(i);
			// cerr << err_i << endl;
			err = max(err, err_i);
		}
		cout << "it: " << it << " err: " << err << endl;
	}

	factor.showmsg();
	return 0;
}


