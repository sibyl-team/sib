// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein
// Author: Anna Paola Muntoni

#ifndef _PARAMS_H
#define _PARAMS_H

typedef long double real_t;

struct Params {
	real_t mu;
	real_t pseed;
	real_t tol;
	int maxit;
	Params(real_t, real_t, real_t, int);
	Params();
};


#endif
