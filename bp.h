// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein
// Author: Alessandro Ingrosso
// Author: Anna Paola Muntoni

#include <vector>
#include <map>
#include <iostream>
#include <omp.h>

#include "params.h"


#ifndef FACTORGRAPH_H
#define FACTORGRAPH_H

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
	Node(int index, Pi const & prob_i, Pr const & prob_g) : index(index), prob_g(prob_g), prob_i(prob_i), f_(0) {}
	int index;
	Pr prob_g;
	Pi prob_i;
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
		std::vector<std::tuple<int, Pi, Pr> > const & individuals = std::vector<std::tuple<int, Pi, Pr> >());
	int find_neighbor(int i, int j) const;
	void add_contact(int i, int j, int t, real_t lambda);
	int add_node(int i);
	void init();
	void set_field(int i, std::vector<int> const & tobs, std::vector<int> const & sobs);
	real_t update(int i, real_t damping);
	void show_graph();
	void show_beliefs(std::ostream &);
	real_t iterate(int maxit, real_t tol, real_t damping);
	real_t iteration(real_t damping);
	real_t loglikelihood() const;
	void show_msg(std::ostream &);
	Params params;
};

std::ostream & operator<<(std::ostream &, FactorGraph const &);

#endif
