#include "params.h"

using namespace std;

Params::Params(shared_ptr<Proba> const & pi,
        shared_ptr<Proba> const & pr,
        real_t pseed,
        real_t psus,
        real_t softconstraint,
        real_t pautoinf) :
		prob_i(pi),
		prob_r(pr),
		pseed(pseed),
		psus(psus),
		softconstraint(softconstraint),
		pautoinf(pautoinf)
	{
		if (pseed + psus > 1)
			throw std::domain_error("pseed and psus are exclusive events but pseed+psus>1");
                if (!pi || !pr)
			throw std::invalid_argument("invalid probability definition");
	}

std::ostream & operator<<(std::ostream & ost, Params const & p)
{
    return ost << "Params("
	<< "prob_i=" << *p.prob_i
        << ",prob_r=" << *p.prob_r
        << ",pseed=" << p.pseed
        << ",psus=" << p.psus
        << ",pautoinf=" << p.pautoinf << ")";
}

std::ostream & operator<<(std::ostream & ost, Proba const & p) { p.print(ost); return ost; }
