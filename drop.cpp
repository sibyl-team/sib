#include "bp.h"
#include <vector>

using namespace std;

real_t set_h(real_t bnew, real_t bold)
{
	if (bnew == 0) {
		if (bold == 0)
			return 1;
		throw;
	}
	return bold / bnew;
}

void drop_sc(FactorGraph & fg, int t, int maxit_sc, real_t tol_sc, real_t damping_bp, real_t damping_sc)
{
	int n = fg.nodes.size();
	vector<vector<real_t>> bts(n);
	vector<vector<real_t>> bgs(n);
	for (int i = 0; i < n; ++i) {
		bts[i] = fg.nodes[i].bt;
		bgs[i] = fg.nodes[i].bg;
	}
	fg.drop_contacts(t);
	for (int it = 0; it < maxit_sc; ++it) {
		real_t err_bp = fg.iteration(damping_bp);
		real_t err_sc = 0.0;
#pragma omp parallel for reduction(max:err_sc)
		for (int i = 0; i < n; ++i) {
			Node & f = fg.nodes[i];
			for (int ti = 0; ti < int(f.bt.size()); ++ti) {
				err_sc = max(err_sc, bts[i][ti] - f.bt[ti]);
				err_sc = max(err_sc, bgs[i][ti] - f.bg[ti]);
				f.ht[ti] *= pow(set_h(f.bt[ti], bts[i][ti]), damping_sc);
				f.hg[ti] *= pow(set_h(f.bg[ti], bgs[i][ti]), damping_sc);
			}
		}
		if (max(err_bp, err_sc) < tol_sc)
			return;
	}
}
