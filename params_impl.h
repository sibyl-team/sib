#include "params.h"

std::ostream & operator<<(std::ostream & ost, Gamma const & g) { return ost << "Gamma(" << g.k << "," << g.mu << ")"; }
std::ostream & operator<<(std::ostream & ost, PriorDiscrete const & d) {
    ost << "PriorDiscrete(";
    for (auto it = d.p.begin(); it < d.p.end() - 1; ++it)
        ost << *it << ",";
    return ost << d.p.back() << ")";
}
std::ostream & operator<<(std::ostream & ost, Exponential const & e) { return ost << "Exp("<< e.mu << ")"; }
std::ostream & operator<<(std::ostream & ost, Uniform const & u) { return ost << "Uniform(" << u.p << ")"; }

template<class Pi, class Pr>
std::ostream & operator<<(std::ostream & ost, Params<Pi,Pr> const & p)
{
    return ost << "Params("
	<< "prob_i=" << p.prob_i
        << ",prob_r=" << p.prob_r
        << ",pseed=" << p.pseed
        << ",psus=" << p.psus << ")";
}

