#include "params.h"
#include <iterator>
#include <algorithm>

using namespace std;

std::ostream & operator<<(std::ostream & ost, Test const & o) {
	return ost << "Test(ps=" << o.ps << ", pi=" << o.pi << ", pr=" << o.pr << ")";
}

Params::Params(shared_ptr<Proba> const & pi,
        shared_ptr<Proba> const & pr,
        real_t pseed,
        real_t psus,
        real_t pautoinf,
	real_t learn_rate) :
		prob_i(pi),
		prob_r(pr),
		pseed(pseed),
		psus(psus),
		pautoinf(pautoinf),
		learn_rate(learn_rate)
{
        if (pseed + psus > 1)
                throw std::domain_error("pseed and psus are exclusive events but pseed+psus>1");
        if (!pi || !pr)
                throw std::invalid_argument("invalid probability definition");
}

PriorDiscrete::PriorDiscrete(Proba const & pr, int T) : Proba(T)
{
        for (int t = 0; t < T; ++t) {
            theta[t] = pr(t);
        }
}

std::ostream & operator<<(std::ostream & ost, Params const & p)
{
    return ost << "Params("
	<< "prob_i=" << *p.prob_i
        << ",prob_r=" << *p.prob_r
        << ",pseed=" << p.pseed
        << ",psus=" << p.psus
        << ",pautoinf=" << p.pautoinf << ")"
        << ",learn_rate=" << p.learn_rate << ")";
}

std::ostream & operator<<(std::ostream & ost, Proba const & p) { p.print(ost); return ost; }


std::ostream & operator<<(std::ostream & ost, RealParams const & p)
{
        ost << "RealParams([";
        for (size_t i = 0; i < p.size(); ++i)
            ost << (i ? ",":"")  << p[i];
        return ost << "])";
}
