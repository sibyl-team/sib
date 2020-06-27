#include "params.h"

using namespace std;

Params::Params(shared_ptr<Proba> const & pi,
        shared_ptr<Proba> const & pr,
        real_t pseed,
        real_t psus,
        real_t fp_rate,
        real_t fn_rate,
        real_t pautoinf,
	real_t mu,
	real_t learn_rate) :
		prob_i(pi),
		prob_r(pr),
		pseed(pseed),
		psus(psus),
		fp_rate(fp_rate),
		fn_rate(fn_rate),
		pautoinf(pautoinf),
		mu(mu),
		learn_rate(learn_rate)
	{
		if (pseed + psus > 1)
			throw std::domain_error("pseed and psus are exclusive events but pseed+psus>1");
                if (!pi || !pr)
			throw std::invalid_argument("invalid probability definition");
	}

PriorDiscrete::PriorDiscrete(Proba const & pr, int T) : p(T)
{
        for (int t = 0; t < T; ++t) {
            p[t] = pr(t);
        }
}

std::ostream & operator<<(std::ostream & ost, Params const & p)
{
    return ost << "Params("
	<< "prob_i=" << *p.prob_i
        << ",prob_r=" << *p.prob_r
        << ",pseed=" << p.pseed
        << ",psus=" << p.psus
        << ",fp_rate=" << p.fp_rate
        << ",fn_rate=" << p.fn_rate
        << ",pautoinf=" << p.pautoinf << ")"
        << ",mu=" << p.mu << ")"
        << ",learn_rate=" << p.learn_rate << ")";
}

std::ostream & operator<<(std::ostream & ost, Proba const & p) { p.print(ost); return ost; }
