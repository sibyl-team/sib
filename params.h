#ifndef PARAMS_H
#define PARAMS_H

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iostream>
#include <exception>
#include <memory>


typedef double real_t;
typedef int times_t;

typedef boost::numeric::ublas::vector<real_t> RealParams;

struct Proba
{
	template<class T>
	Proba(T const & n) : theta(n) {}
	virtual real_t operator()(real_t) const = 0;
	virtual RealParams grad(real_t d) const = 0;
	virtual void print(std::ostream &) const = 0;
	std::istream & operator>>(std::istream & ist) {
		for (int i = 0; i < int(theta.size()); ++i)
			return ist >> theta(i);
		return ist;
	}
	RealParams theta;
};

std::ostream & operator<<(std::ostream & ost, Proba const & p);

struct PriorDiscrete : public Proba
{
	PriorDiscrete(std::vector<real_t> const & p) : Proba(p.size()) { for (int t = 0; t < int(p.size()); ++t) theta(t) = p[t]; }
	PriorDiscrete(Proba const & p, int T);
	real_t operator()(real_t d) const { return d < 0 || d >= int(theta.size()) ? 0.0 : theta(d); }
	RealParams grad(real_t d) const {
		int const T = theta.size();
		RealParams dp(T);
		for (int t = 0; t < T; ++t)
			dp[t] = (-(T-t)*1.0/T + (t <= d));
		return dp;
	}
	void print(std::ostream & ost) const {
		ost << "PriorDiscrete(";
		for (auto it = theta.begin(); it < theta.end() - 1; ++it)
			ost << *it << ",";
		ost << *theta.end() << ")";
	}
};



struct Uniform : public Proba
{
	Uniform(real_t p) : Proba(RealParams(1, p)) {}
	real_t operator()(real_t d) const { return theta(0); }
	RealParams grad(real_t d) const { return RealParams(1, 1.0); }
	void print(std::ostream & ost) const { ost << "Uniform(" << theta(0) << ")"; }
};




struct Exponential : public Proba
{
	Exponential(real_t mu) : Proba(RealParams(1,mu)) {}
	real_t operator()(real_t d) const { return exp(-theta(0)*d); }
	RealParams grad(real_t d) const { return RealParams(1, -d*exp(-theta(0)*d)); }
	void print(std::ostream & ost) const { ost << "Exponential("<< theta(0) << ")"; }
};


struct Gamma : public Proba
{
	Gamma(real_t k, real_t mu) : Proba(2) {}
	real_t operator()(real_t d) const { return boost::math::gamma_q(theta(0),d*theta(1)); }
	RealParams grad(real_t d) const {
  		auto const x = boost::math::differentiation::make_ftuple<real_t, 1, 1>(theta(0), theta(1));
		auto const & xk = std::get<0>(x);
		auto const & xmu = std::get<1>(x);
		auto const f = boost::math::gamma_q(xk, xmu * d);
		RealParams g(2);
		g(0) = f.derivative(1,0);
		g(1) = f.derivative(0,1);
		return g;
	}
	void print(std::ostream & ost) const { ost << "Gamma(" << theta(0) << "," << theta(1) << ")"; }
};


struct Params {
	std::shared_ptr<Proba> prob_i;
	std::shared_ptr<Proba> prob_r;
	real_t pseed;
	real_t psus;
	real_t fp_rate;
	real_t fn_rate;
	real_t pautoinf;
	real_t learn_rate;
	Params(std::shared_ptr<Proba> const & pi, std::shared_ptr<Proba> const & pr, real_t pseed, real_t psus, real_t fp_rate, real_t fn_rate, real_t pautoinf, real_t learn_rate);
};

std::ostream & operator<<(std::ostream &, Params const &);

#endif
