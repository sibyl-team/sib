#include <string>
#include <iostream>

#include "params.h"
#include <getopt.h>

using namespace std;


Params::Params(int & argc, char ** argv) : obs_file("/dev/null"), cont_file("/dev/null"), mu(1.0), pseed(1e-3), tol(1e-3), maxit(100)
{
	int c;
	while ((c = getopt(argc, argv, "s:i:m:o:c:t:h")) != -1 ) {
		switch(c) {
			case 't':
				tol = stod(string(optarg));
				break;
			case 'i':
				maxit = stod(string(optarg));
				break;
			case 'm':
				mu = stod(string(optarg));
				break;
			case 's':
				pseed = stod(string(optarg));
				break;
			case 'o':
				obs_file = optarg;
				break;
			case 'c':
				cont_file = optarg;
				break;
			case 'h':
				cout << "SIR inference, continuous time" << endl;
				cout << "-c : Contact file with format 'i,j,lambdaij,t' " << endl;
				cout << "-o : Observation file with format 'i,state,t' " << endl;
				cout << "-m : mu parameter " << endl;
				cout << "-t : tolerance for convergence " << endl;
				cout << "-i : max iterations for convergence " << endl;
				exit(1);
			default:
				exit(1);
		}
	}
}

