#ifndef PARAMS_H
#define PARAMS_H

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <iostream>
#include <exception>
#include <memory>


using namespace boost::math::differentiation;
typedef double real_t;
typedef int times_t;

struct Proba
{
	virtual real_t operator()(real_t) const = 0;
	virtual void grad_add(real_t d, real_t p) { };
	virtual void grad_init() {};
	virtual void grad_mul(real_t) {};
	virtual void print(std::ostream &) const = 0;
};

std::ostream & operator<<(std::ostream & ost, Proba const & p);

struct PriorDiscrete : public Proba
{
	PriorDiscrete(std::vector<real_t> const & p) : p(p) {}
	PriorDiscrete(Proba const & p, int T);
	real_t operator()(real_t d) const { return d < 0 || d >= int(p.size()) ? 0.0 : p[d]; }
	std::vector<real_t> p;
	std::vector<real_t> dp;
	void grad_add(real_t d, real_t w) {
		int T = p.size();
		for (int t = 0; t < T; ++t)
			dp[t] += (-(T-t)/T + (t <= d)) * w;
	}
	void grad_mul(real_t w) {
		int T = p.size();
		for (int t = 0; t < T; ++t)
			dp[t] *= w;
       	}
	void grad_init() {
		int T = p.size();
		for (int t = 0; t < T; ++t)
			dp[t] = 0;
	}
	void print(std::ostream & ost) const {
	    ost << "PriorDiscrete(";
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
	Exponential(real_t mu) : mu(mu), dmu(0) {}
	real_t mu;
	real_t dmu;
	real_t operator()(real_t d) const { return exp(-mu*d); }
	void grad_add(real_t d, real_t p) { dmu += -d * exp(-mu*d) * p; }
	void grad_mul(real_t p) { dmu *= p; }
	void grad_init() { dmu = 0; }
	std::istream & operator>>(std::istream & ist) { return ist >> mu; }
	void print(std::ostream & ost) const { ost << "Exponential("<< mu << ")"; }
};


struct Gamma : public Proba
{
	real_t k;
	real_t mu;
	real_t dk;
	real_t dmu;
	Gamma(real_t k, real_t mu) : k(k), mu(mu), dk(0), dmu(0) {}
	real_t operator()(real_t d) const { return boost::math::gamma_q(k,d*mu); }
	void grad_init() { dk = 0.0; dmu = 0.0; }
	void grad_add(real_t d, real_t p) {
  		auto const x = make_ftuple<real_t, 1, 1>(k, mu);
		auto const & xk = std::get<0>(x);
		auto const & xmu = std::get<1>(x);
		auto const f = boost::math::gamma_q(xk, xmu * d);
		dk += f.derivative(1,0) * p;
		dmu += f.derivative(0,1) * p;
	}
	void grad_mul(real_t p) { dk *= p; dmu *= p; }
	std::istream & operator>>(std::istream & ist) { return ist >> k >> mu; }
	void print(std::ostream & ost) const { ost << "Gamma(" << k << "," << mu << ")"; }
};


struct Params {
	std::shared_ptr<Proba> prob_i;
	std::shared_ptr<Proba> prob_r;
	real_t pseed;
	real_t psus;
	real_t fp_rate;
	real_t fn_rate;
	real_t pautoinf;
	real_t mu;
	real_t learn_rate;
	Params(std::shared_ptr<Proba> const & pi, std::shared_ptr<Proba> const & pr, real_t pseed, real_t psus, real_t fp_rate, real_t fn_rate, real_t pautoinf, real_t mu, real_t learn_rate);
};

std::ostream & operator<<(std::ostream &, Params const &);

#endif
