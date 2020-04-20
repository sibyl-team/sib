#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <getopt.h>
#include <stdbool.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

#ifndef siblib
#define siblib




class Params {
	public:
		char * obs_file, * cont_file;
		Params () {
			obs_file = 0;
			cont_file = 0;
		}
		int read_input(int & argc, char ** argv) {
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
			return 0;
		}
};


class Factor {
	public: // I don't know what it means...

		struct neigh {
			int index;  // index of the node
			int pos;    // position of the node in neighbors list
			int nij;     // number of contacts (i,j)  + 2 ? useful ?
			vector<int> times; // 
			std::vector<double> msg; // BP msg nij^2 or
			//std::vector< std::vector<double> > msg; // BP msg nij x nij 
		};

		struct node {
			int index;
			std::vector<double> Ti;  // marginals infection times T[ni+2]
			std::vector<double> Gi;  // marginals recovery times	G[ni+2]
			int ni;		   // tot number of contacts
			std::vector<neigh> neighs;	   // list of neighbors
		};

		int N = 0;
		std::vector<node> nodes;
		Params * params;

		Factor(Params * _params): // to be done
			params(_params) {
				build_factor();
			}

		int update_nodes(int i, int j, int t) {

			bool found = false;
			int l;
			for(l = 0; l < int(nodes.size()); l++) {
				if(nodes[l].index == i) {
					found = true;
					break;
				}
			}
			if(found) {
				std::vector<neigh> auxn = nodes[l].neighs;
				bool findj = false;
				int k;
				for (k = 0; k < int(auxn.size()); k++) {
					if(auxn[k].index == j) {
						findj = true;
						break;
					}
				}
				if(findj) {
					nodes[l].neighs[k].times.push_back(t);
					nodes[l].neighs[k].nij++;
				}else {
					neigh auxj;
					auxj.index = j;
					auxj.pos = int(auxn.size());
					auxj.nij = 1;
					auxj.times.push_back(t);
					nodes[l].neighs.push_back(auxj);
					nodes[l].ni++;
				}

			} else {
				N++;
				node aux;
				aux.index = i;
				neigh auxj;
				auxj.index = j;
				auxj.pos = 0;
				auxj.nij = 1;
				auxj.times.push_back(t);
				aux.ni = 1;
				aux.neighs.push_back(auxj);
				nodes.push_back(aux);
			}
			return 0;

		}

		int update_vec_size() {

			for (int i = 0; i < int(nodes.size()); i++) {
				nodes[i].Ti.resize(nodes[i].ni + 2);
				nodes[i].Gi.resize(nodes[i].ni + 2);
				for (int j = 0; j  < int(nodes[i].neighs.size()); j++) {
					int nij = nodes[i].neighs[j].nij;
					nodes[i].neighs[j].msg.resize((nij + 2)*(nij + 2));
				}
			}

			return 0;
		}


		int build_factor() {

			char * obs_file = params->obs_file;
			char * cont_file = params->cont_file;
			string line;

			ifstream obs (obs_file);
			ifstream cont (cont_file);

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
						update_nodes(i,j,t);
						update_nodes(j,i,t);
					}
				}
				cont.close();
			} else {
				fprintf(stderr, "Error opening %s\n", cont_file);
				exit(EXIT_FAILURE);
			}
			update_vec_size();
			fprintf(stderr, "Number of nodes %d\n", N);
			for(int i = 0; i < int(nodes.size()); i++) {
				fprintf(stderr, "### index %d ###\n", nodes[i].index);
				fprintf(stderr, "### in contact with %d nodes\n", nodes[i].ni);
				std::vector<neigh> aux = nodes[i].neighs;
				for (int j = 0; j < int(aux.size()); j++) {
					fprintf(stderr, "# neighbor %d\n", aux[j].index);
					fprintf(stderr, "# in position %d\n", aux[j].pos);
					fprintf(stderr, "# in contact %d times, in t: ", aux[j].nij);
					for(int t = 0; t < int(aux[j].times.size()); t++)
						fprintf(stderr, "%d ", aux[j].times[t]);
					fprintf(stderr, "\n");
				}
			}
			return 0;
		}
};



#endif
