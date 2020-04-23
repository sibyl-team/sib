#ifndef _PARAMS_H
#define _PARAMS_H

typedef long double real_t;

struct Params {
	char const * obs_file;
	char const * cont_file;
	real_t mu;
	real_t pseed;
	real_t tol;
	int maxit;
	Params(int &, char **);
};


#endif
