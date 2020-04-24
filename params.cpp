// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein
// Author: Anna Paola Muntoni

#include <string>
#include <iostream>

#include "params.h"
#include <getopt.h>

using namespace std;


Params::Params() : mu(1.0), pseed(1e-3), tol(1e-3), maxit(100) {}
Params::Params(real_t mu, real_t pseed, real_t tol, int maxit) : mu(mu), pseed(pseed), tol(tol), maxit(maxit) {}

std::ostream & operator<<(std::ostream & o, Params & p)
{
	return o << "Params(mu="
		<< p.mu << ", pseed="
		<< p.pseed << ", tol="
		<< p.tol << ", maxit="
		<< p.maxit << ")" << endl;
}
