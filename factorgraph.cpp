#include <string.h>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cmath>
#include "factorgraph.h"
#include "cavity.h"

using namespace std;

Params::Params(int & argc, char ** argv) : obs_file("/dev/null"), cont_file("/dev/null")
{
	int c;
	while ((c = getopt(argc, argv, "o:c:h")) != -1 ) {
		switch(c) {
			case 'm':
				mu = stod(string(optarg));
				break;
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
				fprintf(stdout, "-m : mu parameter'\n");
				exit(1);
			default:
				exit(1);
		}
	}
}

FactorGraph::FactorGraph(Params const & params) : params(params)
{
	char const * obs_file = params.obs_file;
	char const * cont_file = params.cont_file;

	Tinf = -1;
	string line;

	ifstream obs(obs_file);
	ifstream cont(cont_file);

	int nlines = 0;
	if (cont.is_open()) {
		while (getline(cont,line)) {
			nlines++;
			if (nlines > 1) {
				stringstream s(line);
				int i, j, t;
				char g1, g2, g3;
				real_t lambda;
				s >> i >> g1 >> j >> g2 >> lambda >> g3 >> t;
				cout << i << " " << j << " " << lambda << " " << t << endl;
				add_contact(i, j, t, lambda);
			}
		}
		cont.close();
	} else {
		cerr << "Error opening " << cont_file << endl;
		exit(EXIT_FAILURE);
	}
	finalize();
	showgraph();

}

int FactorGraph::find_neighbor(int i, int j) const
{
	int k = 0;
	for (; k < int(nodes[i].neighs.size()); ++k)
		if (j == nodes[i].neighs[k].index)
			break;
	return k;
}

int FactorGraph::add_node(int i)
{
	map<int,int>::iterator mit = index.find(i);
	if (mit != index.end())
		return mit->second;
	index[i] = nodes.size();
	nodes.push_back(Node(i, params.mu));
	return index[i];
}

void FactorGraph::add_contact(int i, int j, int t, real_t lambda)
{
	Tinf = max(Tinf, t + 1);
	i = add_node(i);
	j = add_node(j);
	int ki = find_neighbor(i, j);
	int kj = find_neighbor(j, i);
	if (ki == int(nodes[i].neighs.size()))
		nodes[i].neighs.push_back(Neigh(j, kj));
	if (kj == int(nodes[j].neighs.size()))
		nodes[j].neighs.push_back(Neigh(i, ki));
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
	for (int k = 0; k < int(F.size()); ++k) {
		if (nodes[i].times.back() != F[k])
			nodes[i].times.push_back(F[k]);
	}
}

void FactorGraph::finalize()
{
	vector<int> F;

	for (int i = 0; i < int(nodes.size()); ++i) {
		finalize_node(i);
		int ntimes = nodes[i].times.size();
		nodes[i].bt.resize(ntimes);
		nodes[i].bg.resize(ntimes);
		nodes[i].ht.resize(ntimes);
		nodes[i].hg.resize(ntimes);
		for (int k = 0; k  < int(nodes[i].neighs.size()); ++k) {
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
			fprintf(stderr, "# neighbor %d\n", nodes[aux[j].index].index);
			fprintf(stderr, "# in position %d\n", aux[j].pos);
			fprintf(stderr, "# in contact %d times, in t: ", int(aux[j].times.size()));
			for (int t = 0; t < int(aux[j].times.size()); t++)
				fprintf(stderr, "%d ", aux[j].times[t]);
			fprintf(stderr, "\n");
		}
	}
}


real_t setmes(vector<real_t> & from, vector<real_t> & to)
{
	int n = from.size();
	real_t s = 0;
	for (int i = 0; i < n; ++i) {
		s += from[i];
	}
	real_t err = 0;
	for (int i = 0; i < n; ++i) {
		from[i] /= s;
		err = max(err, abs(from[i] - to[i]));
		to[i] = from[i];
	}
	return err;
}

int Sij(Node const & f, int j, int sij, int gi) {
	// here gi stands for the ti + gi index
	return f.neighs[j].times[sij] <= f.times[gi] ? sij : f.neighs[j].times.size() - 1;
}



int idx(int sij, int sji, int qj) { return sji + qj * sij; }

real_t prob_obs(Node const & f, int gi, int ti) {
	return exp(-f.mu * (f.times[gi] - f.times[ti])) - exp(-f.mu * (f.times[gi + 1] - f.times[ti]));
}

real_t FactorGraph::update(int i)
{
	int const n = nodes[i].neighs.size();

	Node & f = nodes[i];

	vector<vector<real_t> > UU(n);
	int const qi_ = f.bt.size();

	vector<real_t> ut(qi_);
	vector<real_t> ug(qi_);

	for (int j = 0; j < n; ++j)
		UU[j].resize(f.neighs[j].msg.size());
	// proba tji >= ti for each j
	vector<real_t> C0(n);
	// proba tji > ti for each j
	vector<real_t> C1(n);

	Cavity<real_t> P0(C0, 1., multiplies<real_t>());
	Cavity<real_t> P1(C1, 1., multiplies<real_t>());
	vector<int> min_in(n), min_out(n);
	for (int ti = 0; ti < qi_; ++ti) {
		for (int j = 0; j < n; ++j) {
			Neigh & v = f.neighs[j];
			int const qj = v.times.size();
			min_in[j] = qj - 1;
			min_out[j] = qj - 1;
			for (int s = qj - 1; s >= 0 && v.times[s] >= f.times[ti]; s--) {
				// smallest tji >= ti
				min_in[j] = s;
				if (v.times[s] > f.times[ti]) {
					// smallest tji > ti
					min_out[j] = s;
				}
			}
		}

		for (int gi = ti; gi < qi_; ++gi) {
			fill(C0.begin(), C0.end(), 0.0);
			fill(C1.begin(), C1.end(), 0.0);
			for (int j = 0; j < n; ++j) {
				Neigh & v = f.neighs[j];
				vector<real_t> & h = nodes[v.index].neighs[v.pos].msg;
				int const qj = v.times.size();
				for (int sji = min_in[j]; sji < qj; ++sji) {
					real_t pi = 1;
					for (int s = min_out[j]; s < qj - 1; ++s) {
						int const sij = Sij(f, j, s, gi);
						real_t const p = pi * v.lambdas[s] * h[idx(sij, sji, qj)];
						C0[j] += p;
						if (v.times[sji] > f.times[ti])
							C1[j] += p;
						pi *= 1 - v.lambdas[s];
					}
					int const sij = Sij(f, j, qj - 1, gi);
					real_t const p = pi * h[idx(sij, sji, qj)];
					C0[j] += p;
					if (v.times[sji] > f.times[ti])
						C1[j] += p;
				}
			}

			P0.initialize(C0.begin(), C0.end(), 1.0, multiplies<real_t>());
			P1.initialize(C1.begin(), C1.end(), 1.0, multiplies<real_t>());
			//message to ti
			// FIX HERE
			real_t g_prob = prob_obs(f, gi, ti);
			real_t a = g_prob  * (ti == 0 || ti == qi_ - 1 ? P0.full() : P0.full() - P1.full());
			ug[gi] += f.bt[ti] * a;
			ut[ti] += f.bg[gi] * a;

			//messages to sij, sji
			for (int j = 0; j < n; ++j) {
				Neigh & v = f.neighs[j];
				int const qj = v.times.size();
				real_t const p0 = P0[j];
				real_t const p01 = p0 - P1[j];
				for (int sji = min_in[j]; sji < qj; ++sji) {
					real_t pi = f.ht[ti] * f.hg[gi] * g_prob * (ti == 0 || v.times[sji] == f.times[ti] ? p0 : p01);
					for (int s = min_out[j]; s < qj - 1; ++s) {
						real_t & Uij = UU[j][idx(Sij(f, j, s, gi), sji, qj)];
						Uij += pi * v.lambdas[s];
						pi *= 1 - v.lambdas[s];
					}
					UU[j][idx(Sij(f, j, qj - 1, gi), sji, qj)] += pi;
				}
			}
		}
	}

	real_t diff = max(setmes(ut, f.bt), setmes(ug, f.bg));
	for (int j = 0; j < n; ++j)
		diff = max(diff, setmes(UU[j], f.neighs[j].msg));

	return diff;

}

