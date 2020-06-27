// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein
// Author: Anna Paola Muntoni

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
#include "bp.h"

using namespace std;

pair<vector<tuple<int,int,times_t,real_t> >, vector<tuple<int,int,times_t> > >
read_files(char const * cont_file, char const * obs_file)
{
	string line;

	ifstream cont(cont_file);
	if (!cont.is_open()) {
		cerr << "Error opening " << cont_file << endl;
		exit(EXIT_FAILURE);
	}

	auto contacts = vector<tuple<int,int,times_t,real_t> >();
	auto observations = vector<tuple<int,int,times_t> >();

	int nlines = 0;
	while (getline(cont, line)) {
		nlines++;
		if (nlines > 1) { // waring: firrst line (header) is silently ignored!!
			stringstream s(line);
			int i, j, t;
			char g1, g2, g3;
			real_t lambda;
			s >> i >> g1 >> j >> g2 >> t >> g3 >> lambda;
			//cout << i << " " << j << " " << lambda << " " << t << endl;
			contacts.push_back(make_tuple(i,j,t,lambda));
		}
	}
	cont.close();

	ifstream obs(obs_file);
	if (!obs.is_open()) {
		cerr << "Error opening " << obs_file << endl;
		exit(EXIT_FAILURE);
	}
	nlines = 0;
	while (getline(obs,line)) {
		nlines++;
		if(nlines > 1) {
			stringstream s(line);
			int i, state, t;
			char g1, g2;
			s >> i >> g1 >> state >> g2 >> t;
			//cout << i << state << t << endl;
			observations.push_back(make_tuple(i,state,t));
		}
	}
	return make_pair(contacts,observations);
}

int main(int argc, char ** argv)
{
	char const * obs_file = "/dev/null";
	char const * cont_file = "/dev/null";
	int c;

	real_t tol = 1e-3;
	real_t mu = 0.5;
	real_t pseed = 0.1;
	int maxit = 100;

	while ((c = getopt(argc, argv, "j:s:i:m:o:c:t:hv")) != -1 ) {
		switch(c) {
			case 'j':
				omp_set_num_threads(stod(optarg));
				break;
			case 't':
				tol = stod(optarg);
				break;
			case 'i':
				maxit = stoi(optarg);
				break;
			case 'm':
				mu = stod(optarg);
				break;
			case 's':
				pseed = stoi(optarg);
				break;
			case 'o':
				obs_file = optarg;
				break;
			case 'c':
				cont_file = optarg;
				break;
			case 'v':
				cout << "SIBYL version " << VERSION << endl;
				return 0;
			case 'h':
				cout << "SIR inference, continuous time" << endl;
				cout << "-c : Contact file with format 'i,j,t,lambdaij' " << endl;
				cout << "-o : Observation file with format 'i,state,t' " << endl;
				cout << "-m : mu parameter " << endl;
				cout << "-t : tolerance for convergence " << endl;
				cout << "-i : max iterations for convergence " << endl;
				exit(1);
			default:
				exit(1);
		}
	}

	auto prob_i = shared_ptr<Proba>(new Uniform(1.0));
	auto prob_r = shared_ptr<Proba>(new Exponential(mu));
	Params p(prob_i, prob_r, pseed, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0);
	auto co = read_files(cont_file, obs_file);
	FactorGraph factor(p, get<0>(co), get<1>(co));
	factor.iterate(maxit, tol, 0.0);
	factor.show_beliefs(cout);
	return 0;
}


