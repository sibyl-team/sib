#include "params.h"

std::ostream & operator<<(std::ostream & ost, Gamma const & g) { return ost << "Gamma(" << g.k << "," << g.mu << ")"; }
std::ostream & operator<<(std::ostream & ost, Exponential const & e) { return ost << "Exp("<< e.mu << ")"; }
std::ostream & operator<<(std::ostream & ost, Uniform const & u) { return ost << "Uniform(" << u.p << ")"; }

std::ostream & operator<<(std::ostream & ost, Params const & p)
{
    return ost << "Params("
	<< "prob_i=" << p.prob_i
        << ",prob_r=" << p.prob_r
        << ",pseed=" << p.pseed
        << ",psus=" << p.psus << ")";
}

