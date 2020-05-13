// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein
// Author: Alessandro Ingrosso
// Author: Anna Paola Muntoni

#include <vector>
#include <iostream>
#include <memory>
#include <omp.h>

#include "params.h"


#ifndef FACTORGRAPH_H
#define FACTORGRAPH_H

struct Mes : public std::vector<real_t>
{
	Mes() : qj(0) {}
	Mes(size_t qj) : vector<real_t>(qj*qj), qj(qj) {}
	void clear() { std::fill(begin(), end(), 0.0); }
	size_t dim() const { return qj;}
	real_t & operator()(int sji, int sij) { return operator[](qj * sij + sji); }
	real_t operator()(int sji, int sij) const { return operator[](qj * sij + sji); }
	size_t qj;
};

struct Neigh {
	Neigh(int index, int pos) : index(index), pos(pos) { omp_init_lock(&lock_); }
	int index;  // index of the node
	int pos;    // position of the node in neighbors list
	std::vector<int> t; // time index of contacts
	std::vector<real_t> lambdas; // transmission probability
	Mes msg; // BP msg nij^2 or
	void lock() const { omp_set_lock(&lock_); }
	void unlock() const { omp_unset_lock(&lock_); }
	mutable omp_lock_t lock_;
};

struct Node {
	Node(std::shared_ptr<Proba> prob_i, std::shared_ptr<Proba> prob_r) :
		prob_i(prob_i),
		prob_r(prob_r),
		f_(0)
	{}
	std::shared_ptr<Proba> prob_i;
	std::shared_ptr<Proba> prob_r;
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
	FactorGraph(Params const & params,
		std::vector<std::tuple<int,int,int,real_t> > const & contacts,
		std::vector<std::tuple<int, int, int> > const & obs,
		std::vector<std::tuple<int, std::shared_ptr<Proba>, std::shared_ptr<Proba>> > const & individuals 
			= std::vector<std::tuple<int, std::shared_ptr<Proba>, std::shared_ptr<Proba>>>());
	int find_neighbor(int i, int j) const;
	void add_contact(int i, int j, int t, real_t lambda);
	int get_node(int i);
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
