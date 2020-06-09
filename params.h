#ifndef PARAMS_H
#define PARAMS_H

#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <exception>
#include <memory>


typedef double real_t;
typedef int times_t;

struct Proba
{
	virtual real_t operator()(real_t) const = 0;
	virtual real_t operator()(real_t d, real_t lambda) const { return operator()(d)*lambda; }
	virtual void print(std::ostream &) const = 0;
};

std::ostream & operator<<(std::ostream & ost, Proba const & p);

struct PriorDiscrete : public Proba
{
	PriorDiscrete(std::vector<real_t> const & p) : p(p) {}
	real_t operator()(real_t d) const { return d < 0 || d >= int(p.size()) ? 0.0 : p[d]; }
	std::vector<real_t> p;
	void print(std::ostream & ost) const {
	    ost << "PriorDiscrete(";
	    for (auto it = p.begin(); it < p.end() - 1; ++it)
		ost << *it << ",";
	    ost << p.back() << ")";
	}
};


struct ExpDiscrete : public Proba
{
	ExpDiscrete(std::vector<real_t> const & p) : p(p) {}
	real_t operator()(real_t d) const { return d < 0 || d >= int(p.size()) ? 0.0 : p[d]; }
	real_t operator()(real_t d, real_t lambda) const { return 1.0-std::exp(-operator()(d)*lambda); }
	std::vector<real_t> p;
	void print(std::ostream & ost) const {
	    ost << "ExpDiscrete(";
	    for (auto it = p.begin(); it < p.end() - 1; ++it)
		ost << *it << ",";
	    ost << p.back() << ")";
	}
};


struct Uniform : public Proba
{
	Uniform(real_t p) : p(p) {}
	real_t p;
	real_t operator()(real_t d) const { return p; }
	std::istream & operator>>(std::istream & ist) { return ist >> p; }
	void print(std::ostream & ost) const { ost << "Uniform(" << p << ")"; }
};




struct Exponential : public Proba
{
	Exponential(real_t mu) : mu(mu) {}
	real_t mu;
	real_t operator()(real_t d) const { return exp(-mu*d); }
	std::istream & operator>>(std::istream & ist) { return ist >> mu; }
	void print(std::ostream & ost) const { ost << "Exponential("<< mu << ")"; }
};


struct Gamma : public Proba
{
	real_t k;
	real_t mu;
	Gamma(real_t k, real_t mu) : k(k), mu(mu) {}
	real_t operator()(real_t d) const { return boost::math::gamma_q(k,d*mu); }
	std::istream & operator>>(std::istream & ist) { return ist >> k >> mu; }
	void print(std::ostream & ost) const { ost << "Gamma(" << k << "," << mu << ")"; }
};


struct Params {
	std::shared_ptr<Proba> prob_i;
	std::shared_ptr<Proba> prob_r;
	real_t pseed;
	real_t psus;
	real_t softconstraint;
	real_t pautoinf;
	Params(std::shared_ptr<Proba> const & pi, std::shared_ptr<Proba> const & pr, real_t pseed, real_t psus, real_t softconstraint, real_t pautoinf) :
		prob_i(pi),
		prob_r(pr),
		pseed(pseed),
		psus(psus),
		softconstraint(softconstraint),
		pautoinf(pautoinf)
	{
		if (pseed + psus > 1)
			throw std::domain_error("pseed and psus are exclusive events but pseed+psus>1");
	}
};

std::ostream & operator<<(std::ostream &, Params const &);

#endif
