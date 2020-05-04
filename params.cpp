#include "params.h"


std::ostream & operator<<(std::ostream & ost, Params const & p)
{
    return ost << "Params("
	<< "prob_i=" << *p.prob_i
        << ",prob_r=" << *p.prob_r
        << ",pseed=" << p.pseed
        << ",psus=" << p.psus << ")";
}

std::ostream & operator<<(std::ostream & ost, Proba const & p) { p.print(ost); return ost; }
