#ifndef PARAMS_H
#define PARAMS_H

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/differentiation/autodiff.hpp>
#include <valarray>

#include <iostream>
#include <exception>
#include <memory>


typedef double real_t;
typedef int times_t;

typedef std::valarray<real_t> RealParams;

struct Proba
{
	template<class T>
	Proba(T const & n) : theta(n) {}
	virtual real_t operator()(real_t) const = 0;
	virtual void grad(RealParams & dtheta, real_t d) const = 0;
	virtual void print(std::ostream &) const = 0;
	std::istream & operator>>(std::istream & ist) {
		for (int i = 0; i < int(theta.size()); ++i)
			return ist >> theta[i];
		return ist;
	}
	RealParams theta;
};

std::ostream & operator<<(std::ostream & ost, Proba const & p);
std::ostream & operator<<(std::ostream & ost, RealParams const & p);

struct PriorDiscrete : public Proba
{
	// PriorDiscrete(std::vector<real_t> const & p) : Proba(p.size()) { for (size_t t = 0; t < p.size(); ++t) theta[t] = p[t]; }
	PriorDiscrete(RealParams const & p) : Proba(p) {}
	PriorDiscrete(Proba const & p, int T);
	real_t operator()(real_t d) const { return d < 0 || d >= int(theta.size()) ? 0.0 : theta[d]; }
	void grad(RealParams & dtheta, real_t d) const {
		for (auto & x : dtheta)
			x = 0.0;
		if (d < dtheta.size())
			dtheta[d] = 1.0;
	}
	void print(std::ostream & ost) const {
		ost << "PriorDiscrete(" << theta << ")";
	}
};

struct PiecewiseLinear : public Proba
{
	PiecewiseLinear(RealParams const & p, real_t step) : Proba(p), step(step) {}
	real_t operator()(real_t d) const {
		real_t const x = d / step;
		if (x < 0 || x > theta.size() - 1)
			return 0;
		int const k = x;
		if (k == x)
			return theta[k];
		return (k + 1 - x) * theta[k] + (x - k) * theta[k + 1];
	}
	void grad(RealParams & dtheta, real_t d) const {
		std::fill(std::begin(dtheta), std::end(dtheta), real_t(0.0));
		real_t const x = d / step;
		if (x < 0 || x > theta.size() - 1)
			return;
		int const k = x;
		if (k == x) {
			dtheta[k] = 1.0;
		} else {
			dtheta[k] = k + 1 - x;
			dtheta[k + 1] = x - k;
		}
	}
	void print(std::ostream & ost) const {
		ost << "PieceWiseLinear(" << theta << ")";
	}
	real_t step;
};



struct Cached : public Proba
{
	Cached(std::shared_ptr<Proba> const & prob, int T) : Proba(prob->theta), prob(prob), p(T), zero(0.0, prob->theta.size()), dp(T, zero)  {
		recompute();
	}
	std::shared_ptr<Proba> prob;
	std::vector<real_t> p;
	RealParams const zero;
	std::vector<RealParams> dp;
	real_t operator()(real_t d) const { return d < 0 || d >= int(p.size()) ? 0.0: p[d]; }
	void grad(RealParams & dtheta, real_t d) const { dtheta = d < 0 || d >= int(dp.size()) ? zero : dp[d]; }
	void recompute() {
		prob->theta = theta;
		for (size_t d = 0; d < p.size(); ++d) {
			p[d] = (*prob)(d);
			prob->grad(dp[d], d);
		}
	}
	void print(std::ostream & ost) const { ost << "Cached(" << prob << ")"; }
};

struct Uniform : public Proba
{
	Uniform(real_t p) : Proba(RealParams({p})) {}
	real_t operator()(real_t d) const { return theta[0]; }
	void grad(RealParams & dtheta, real_t d) const { dtheta[0] = 1.0; }
	void print(std::ostream & ost) const { ost << "Uniform(" << theta[0] << ")"; }
};




struct Exponential : public Proba
{
	Exponential(real_t mu) : Proba(RealParams({mu})) {}
	real_t operator()(real_t d) const { return exp(-theta[0]*d); }
	void grad(RealParams & dtheta, real_t d) const { dtheta[0]= -d*exp(-theta[0]*d); }
	void print(std::ostream & ost) const { ost << "Exponential("<< theta[0] << ")"; }
};


struct Gamma : public Proba
{
	Gamma(real_t k, real_t mu) : Proba(RealParams({k,mu})) {}
	real_t operator()(real_t d) const { return boost::math::gamma_q(theta[0], d * theta[1]); }
	void grad(RealParams & dtheta, real_t d) const {
		if (!d) {
			dtheta[0] = 0.0;
			dtheta[1] = 0.0;
			return;
		}
  		auto const x = boost::math::differentiation::make_ftuple<real_t, 1, 1>(theta[0], theta[1]);
		auto const & xk = std::get<0>(x);
		auto const & xmu = std::get<1>(x);
		auto const f = boost::math::gamma_q(xk, xmu * d);
		dtheta[0] = f.derivative(1,0);
		dtheta[1] = f.derivative(0,1);
	}
	void print(std::ostream & ost) const { ost << "Gamma(" << theta[0] << "," << theta[1] << ")"; }
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
