// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein
// Author: Alessandro Ingrosso
// Author: Anna Paola Muntoni

#include <vector>
#include <map>
#include <iostream>
#include "omp.h"

#ifndef FACTORGRAPH_H
#define FACTORGRAPH_H

typedef long double real_t;

struct Params {
	real_t mu;
	real_t k;
	real_t pseed;
	real_t damping;
	// for k=1.0, the gamma becomes exponential
	Params() : mu(0.5), k(1.0), pseed(1e-3), damping(0.) {}
	Params(real_t mu, real_t pseed, real_t damping) : mu(mu), k(1.0), pseed(pseed), damping(damping) {}
	Params(real_t mu, real_t k, real_t pseed, real_t damping) : mu(mu), k(k), pseed(pseed), damping(damping) {}
};

std::ostream & operator<<(std::ostream &, Params const &);

struct Neigh {
	Neigh(int index, int pos) : index(index), pos(pos) {}
	int index;  // index of the node
	int pos;    // position of the node in neighbors list
	std::vector<int> times; // times of contacts
	std::vector<real_t> lambdas; // times of contacts
	std::vector<real_t> msg; // BP msg nij^2 or
	omp_lock_t lock_;
};

struct Node {
	Node(int index, real_t k, real_t mu) : index(index), k_(k), mu_(mu), f_(0) {}
	int index;
	real_t k_;
	real_t mu_;
	real_t prob_g(real_t delta) const;
	std::vector<int> tobs;
	std::vector<int> obs;
	std::vector<int> times;
	std::vector<real_t> bt;  // marginals infection times T[ni+2]
	std::vector<real_t> bg;  // marginals recovery times G[ni+2]
	std::vector<real_t> ht;  // message infection times T[ni+2]
	std::vector<real_t> hg;  // message recovery times G[ni+2]
	std::vector<Neigh> neighs;	   // list of neighbors
	real_t f_;
};


class FactorGraph {
public:
	int Tinf;
	std::vector<Node> nodes;
	std::map<int, int> index;
	FactorGraph(Params const & params,
		std::vector<std::tuple<int,int,int,real_t> > const & contacts,
		std::vector<std::tuple<int, int, int> > const & obs,
		std::vector<std::tuple<int, real_t, real_t> > const & individuals = std::vector<std::tuple<int, real_t, real_t> >());
	int find_neighbor(int i, int j) const;
	void add_contact(int i, int j, int t, real_t lambda);
	int add_node(int i);
	void add_obs(int i, int state, int t);
	void init();
	void set_field(int i);
	real_t update(int i);
	void show_graph();
	void show_beliefs(std::ostream &);
	real_t iterate(int maxit, real_t tol);
	real_t iteration();
	real_t loglikelihood() const;
	void show_msg(std::ostream &);

	std::map<int, std::vector<real_t> > get_tbeliefs();
	std::map<int, std::vector<real_t> > get_gbeliefs();
	Params params;
};


#endif
