#include <vector>
#include <map>

#ifndef FACTORGRAPH_H
#define FACTORGRAPH_H



struct Params {
	char const * obs_file;
	char const * cont_file;
	Params(int &, char **);
};

struct Neigh {
	Neigh(int index, int pos) : index(index), pos(pos) {}
	int index;  // index of the node
	int pos;    // position of the node in neighbors list
	std::vector<int> times; // times of contacts
	std::vector<float> lambdas; // times of contacts
	std::vector<double> msg; // BP msg nij^2 or
};

struct Node {
	Node(int index) : index(index) {}
	int index;
	std::vector<int> times;
	std::vector<double> Ti;  // marginals infection times T[ni+2]
	std::vector<double> Gi;  // marginals recovery times G[ni+2]
	std::vector<Neigh> neighs;	   // list of neighbors
};


class FactorGraph {
public:
	int Tinf;
	std::vector<Node> nodes;
	std::map<int, int> idx;
	FactorGraph(Params const & params);
	int find_neighbor(int i, int j) const;
	void add_contact(int i, int j, int t, float lambda);
	int add_node(int i);
	void finalize_node(int i);
	void finalize();
	void showgraph();
};



#endif
