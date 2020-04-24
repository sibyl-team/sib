// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein
// Author: Anna Paola Muntoni

#include <string>
#include <iostream>

#include "params.h"
#include <getopt.h>

using namespace std;


Params::Params() : mu(1.0), pseed(1e-3) {}
Params::Params(real_t mu, real_t pseed) : mu(mu), pseed(pseed) {}

std::ostream & operator<<(std::ostream & o, Params & p)
{
	return o << "Params(mu="
		<< p.mu << ", pseed="
		<< p.pseed << ")" << endl;
}
