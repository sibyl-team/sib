#include "factorgraph.h"
#include <algorithm>

using namespace std;

Params::Params(int & argc, char ** argv) {
	int c;
	while ((c = getopt(argc, argv, "o:c:h")) != -1 ) {
		switch(c) {
			case 'o':
				obs_file = optarg;
				break;
			case 'c':
				cont_file = optarg;
				break;
			case 'h':
				fprintf(stdout, "SIR inference, continuous time\n");
				fprintf(stdout, "-c : Contact file with format 'i,j,lambdaij,t'\n");
				fprintf(stdout, "-o : Observation file with format 'i,state,t'\n");
				exit(1);
			default:
				exit(1);
		}
	}
}

FactorGraph::FactorGraph(Params const & params)
{
	char const * obs_file = params.obs_file;
	char const * cont_file = params.cont_file;
	string line;

	ifstream obs(obs_file);
	ifstream cont(cont_file);

	int nlines = 0;
	if (cont.is_open()) {
		while (getline(cont,line)) {
			nlines++;
			if(nlines > 1) {
				stringstream s(line);
				int i, j, t;
				char g1, g2, g3;
				double lambda;
				s >> i >> g1 >> j >> g2 >> lambda >> g3 >> t;
				fprintf(stdout, "%d %d %f %d\n", i,j,lambda,t);
				add_contact(i, j, t, lambda);
			}
		}
		cont.close();
	} else {
		fprintf(stderr, "Error opening %s\n", cont_file);
		exit(EXIT_FAILURE);
	}
	finalize();
	showgraph();

}

int FactorGraph::find_neighbor(int i, int j) const
{
	int k = 0;
	for (; k < int(nodes[i].neighs.size()); ++j)
		if (j == nodes[i].neighs[k].index)
			break;
	return k;
}

int FactorGraph::add_node(int i)
{
	map<int,int>::iterator mit = idx.find(i);
	if (mit != idx.end())
		return mit->second;
	idx[i] = nodes.size();
	nodes.push_back(Node(i));
	return idx[i];
}

void FactorGraph::add_contact(int i, int j, int t, float lambda)
{
	i = add_node(i);
	j = add_node(j);
	int ki = find_neighbor(i, j);
	int kj = find_neighbor(j, i);
	if (ki == int(nodes[i].neighs.size()))
		nodes[i].neighs.push_back(Neigh(i, kj));
	if (kj == int(nodes[j].neighs.size()))
		nodes[j].neighs.push_back(Neigh(j, ki));
	nodes[i].neighs[ki].times.push_back(t);
	nodes[i].neighs[ki].lambdas.push_back(lambda);
	nodes[j].neighs[kj].times.push_back(t);
	nodes[j].neighs[kj].lambdas.push_back(lambda);
}

void FactorGraph::finalize_node(int i)
{
	vector<int> F;
	for (int k = 0; k < int(nodes[i].neighs.size()); ++k) {
		vector<int> const & tij = nodes[i].neighs[k].times;
		F.insert(F.end(), tij.begin(), tij.end());
	}
	sort(F.begin(), F.end());
	F.push_back(Tinf);
	F.push_back(numeric_limits<int>::max());
	nodes[i].times.push_back(-1);
	for (int i = 0; i < int(F.size()); ++i) {
		if (nodes[i].times.back() != F[i])
			nodes[i].times.push_back(F[i]);
	}
}

void FactorGraph::finalize()
{
	vector<int> F;

	for (int i = 0; i < int(nodes.size()); i++) {
		int ntimes = nodes[i].times.size();
		nodes[i].Ti.resize(ntimes);
		nodes[i].Gi.resize(ntimes);
		for (int k = 0; k  < int(nodes[i].neighs.size()); k++) {
			int nij = nodes[i].neighs[k].times.size();
			nodes[i].neighs[k].msg.resize((nij + 2)*(nij + 2));
		}
	}
}


void FactorGraph::showgraph()
{
	fprintf(stderr, "Number of nodes %d\n", int(nodes.size()));
	for(int i = 0; i < int(nodes.size()); i++) {
		fprintf(stderr, "### index %d ###\n", nodes[i].index);
		fprintf(stderr, "### in contact with %d nodes\n", int(nodes[i].neighs.size()));
		vector<Neigh> const & aux = nodes[i].neighs;
		for (int j = 0; j < int(aux.size()); j++) {
			fprintf(stderr, "# neighbor %d\n", aux[j].index);
			fprintf(stderr, "# in position %d\n", aux[j].pos);
			fprintf(stderr, "# in contact %d times, in t: ", int(aux[j].times.size()));
			for (int t = 0; t < int(aux[j].times.size()); t++)
				fprintf(stderr, "%d ", aux[j].times[t]);
			fprintf(stderr, "\n");
		}
	}
}
