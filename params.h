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

struct Test
{
	Test(real_t ps, real_t pi, real_t pr) : ps(ps), pi(pi), pr(pr) {}
	real_t ps;
	real_t pi;
	real_t pr;
};


std::ostream & operator<<(std::ostream & ost, Test const & o);

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
	virtual void set_theta(RealParams const & newtheta) { theta = newtheta; }
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

	PiecewiseLinear(Proba const & p, int T, real_t step = 1.0) : Proba(T), step(step) {
		for (int t = 0; t < T; ++t)
			theta[t] = p(t * step);
	}

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
		update();
	}
	std::shared_ptr<Proba> prob;
	std::vector<real_t> p;
	RealParams const zero;
	std::vector<RealParams> dp;
	real_t operator()(real_t d) const { return d < 0 || d >= int(p.size()) ? 0.0: p[d]; }
	void grad(RealParams & dtheta, real_t d) const { dtheta = d < 0 || d >= int(dp.size()) ? zero : dp[d]; }
	void update() {
		prob->set_theta(theta);
		for (size_t d = 0; d < p.size(); ++d) {
			p[d] = (*prob)(d);
			prob->grad(dp[d], d);
		}
	}

	virtual RealParams get_theta() const { return theta; }
	void set_theta(RealParams const & newtheta) {
		theta = newtheta;
		update();
	}
	void print(std::ostream & ost) const { ost << "Cached(" << *prob << ",T=" << p.size() << ")"; }
};

struct Scaled : public Proba
{
	Scaled(std::shared_ptr<Proba> const & prob, real_t scale) : Proba(prob->theta.size() + 1), prob(prob) {
		for (size_t i = 0; i < theta.size() - 1; ++i)
			theta[i] = prob->theta[i];
		theta[theta.size() - 1] = scale;
	}
	real_t operator()(real_t d) const { return prob->operator()(d) * theta[theta.size() - 1]; }
	void grad(RealParams & dtheta, real_t d) const {
		prob->grad(dtheta, d);
		dtheta *= theta[theta.size() - 1];
		dtheta[dtheta.size() - 1] = prob->operator()(d);
	}

	void set_theta(RealParams const & newtheta) {
		theta = newtheta;
		prob->set_theta(RealParams(&theta[0], prob->theta.size()));
	}
	void print(std::ostream & ost) const { ost << "Scaled(" << *prob << ",scale=" << theta[theta.size()-1] << ")"; }
	std::shared_ptr<Proba> prob;
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

struct UnnormalizedGammaPDF : public Proba
{
	UnnormalizedGammaPDF(real_t k, real_t mu) : Proba(RealParams({k,mu})) {}
	real_t operator()(real_t d) const { return d ? exp(-theta[1] * d + (theta[0]-1) * log(d)) : 0.0; }
	void grad(RealParams & dtheta, real_t d) const {
		real_t const p = operator()(d);
		dtheta[0] = d ? log(d) * p : 0.0;
		dtheta[1] = -d * p;
	}
	void print(std::ostream & ost) const { ost << "UnnormalizedGammaPDF(" << theta[0] << "," << theta[1] << ")"; }
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
		auto const f = boost::math::gamma_q(std::get<0>(x), std::get<1>(x) * d);
		dtheta[0] = f.derivative(1,0);
		dtheta[1] = f.derivative(0,1);
	}
	void print(std::ostream & ost) const { ost << "Gamma(" << theta[0] << "," << theta[1] << ")"; }
};

struct ConstantRate : public Proba
{
	ConstantRate(real_t gamma, real_t Dt) : Proba(RealParams({gamma})), Dt(Dt) {}
	real_t operator()(real_t d) const { return -expm1(-theta[0]*Dt); }
	void grad(RealParams & dtheta, real_t d) const { dtheta[0] = Dt*exp(-theta[0]*Dt); }
	void print(std::ostream & ost) const { ost << "ConstatRate(" << theta[0] << "," << Dt << ")"; }
	real_t Dt;
};


struct PDF : public Proba
{
	PDF(std::shared_ptr<Proba> const & prob) : Proba(prob->theta), prob(prob) {}
	real_t operator()(real_t d) const { return prob->operator()(d) - prob->operator()(d + 1); }
	void grad(RealParams & dtheta, real_t d) const {
		prob->grad(dtheta, d);
		RealParams dtheta1(dtheta.size());
		prob->grad(dtheta1, d + 1);
		dtheta -= dtheta1;
	}

	void set_theta(RealParams const & newtheta) {
		theta = newtheta;
		prob->set_theta(theta);
	}
	void print(std::ostream & ost) const { ost << "PDF(" << *prob << ")"; }
	std::shared_ptr<Proba> prob;
};


struct Params {
	std::shared_ptr<Proba> prob_i;
	std::shared_ptr<Proba> prob_r;
	std::vector<std::shared_ptr<Test>> obs;
	std::shared_ptr<Test> fakeobs;
	real_t pseed;
	real_t psus;
	real_t pautoinf;
	real_t learn_rate;
	Params(std::shared_ptr<Proba> const & pi, std::shared_ptr<Proba> const & pr, real_t pseed, real_t psus, real_t fp_rate, real_t fn_rate, real_t pautoinf, real_t learn_rate);
};

std::ostream & operator<<(std::ostream &, Params const &);

#endif
