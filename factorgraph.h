#include <vector>
#include <map>

#ifndef FACTORGRAPH_H
#define FACTORGRAPH_H


typedef long double real_t;

struct Params {
	char const * obs_file;
	char const * cont_file;
	real_t mu;
	Params(int &, char **);
};

struct Neigh {
	Neigh(int index, int pos) : index(index), pos(pos) {}
	int index;  // index of the node
	int pos;    // position of the node in neighbors list
	std::vector<int> times; // times of contacts
	std::vector<real_t> lambdas; // times of contacts
	std::vector<real_t> msg; // BP msg nij^2 or
};

struct Node {
	Node(int index, real_t mu) : index(index), mu(mu) {}
	int index;
	real_t mu;
	std::vector<int> times;
	std::vector<real_t> bt;  // marginals infection times T[ni+2]
	std::vector<real_t> bg;  // marginals recovery times G[ni+2]
	std::vector<real_t> ht;  // message infection times T[ni+2]
	std::vector<real_t> hg;  // message recovery times G[ni+2]
	std::vector<Neigh> neighs;	   // list of neighbors
};


class FactorGraph {
public:
	int Tinf;
	std::vector<Node> nodes;
	std::map<int, int> index;
	FactorGraph(Params const & params);
	int find_neighbor(int i, int j) const;
	void add_contact(int i, int j, int t, real_t lambda);
	int add_node(int i);
	void finalize_node(int i);
	void finalize();
	real_t update(int i);
	void showgraph();
	Params params;
};



#endif
