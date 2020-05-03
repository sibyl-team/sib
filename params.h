#ifndef PARAMS_H
#define PARAMS_H

#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <exception>


typedef double real_t;

class Uniform;
class Exponential;
class Gamma;
class ExpGammaInc;

typedef Gamma Pr;
typedef ExpGammaInc Pi;

struct Uniform
{
	Uniform(real_t p) : p(p) {}
	real_t p;
	real_t operator()(real_t d) const { return p; }
	std::istream & operator>>(std::istream & ist) { return ist >> p; }
};


std::ostream & operator<<(std::ostream & ost, Uniform const & u);


struct Exponential
{
	Exponential(real_t mu) : mu(mu) {}
	real_t mu;
	real_t operator()(real_t d) const { return exp(-mu*d); }
	std::istream & operator>>(std::istream & ist) { return ist >> mu; }
};

std::ostream & operator<<(std::ostream & ost, Exponential const & e);


struct Gamma
{
	real_t k;
	real_t mu;
	Gamma(real_t k, real_t mu) : k(k), mu(mu) {}
	real_t operator()(real_t d) const { return 1-boost::math::gamma_p(k,d*mu); }
	std::istream & operator>>(std::istream & ist) { return ist >> k >> mu; }
};

std::ostream & operator<<(std::ostream & ost, Gamma const & g);


struct ExpGammaInc
{
	real_t k;
	real_t mu;
	real_t probs[200];
	ExpGammaInc(real_t k, real_t mu) : k(k), mu(mu) {
		real_t sum = 0;
		for (int d = 0; d < 200; d++) {
			double l = boost::math::gamma_p(k, d * mu);
			if (d > 0)
				l -= boost::math::gamma_p(k, (d - 1) * mu);
			real_t p = 1. - exp(-l);
			probs[d] = p;
			sum += p;
		}
		for (int d = 0; d < 200; d++) {
			probs[d] /= sum;
		}
	}
	real_t operator()(real_t d) const { return probs[int(d)]; }
	std::istream & operator>>(std::istream & ist) { return ist >> k >> mu; }
};

std::ostream & operator<<(std::ostream & ost, ExpGammaInc const & g);


struct Params {
	Pi prob_i;
	Pr prob_r;
	real_t pseed;
	real_t psus;
	Params(Pi const & pi, Pr const & pr, real_t pseed, real_t psus) : prob_i(pi), prob_r(pr), pseed(pseed), psus(psus) {
		if (pseed + psus > 1)
			throw std::domain_error("pseed and psus are exclusive events but pseed+psus>1");
	}
};

std::ostream & operator<<(std::ostream &, Params const &);

#endif
