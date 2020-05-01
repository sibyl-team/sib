// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein
// Author: Alessandro Ingrosso
// Author: Anna Paola Muntoni
// Author: Indaco Biazzo



#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <functional>
#include <algorithm>
#include <cassert>
#include <tuple>
#include <exception>
#include "bp.h"
#include "cavity.h"

using namespace std;




FactorGraph::FactorGraph(Params const & params,
		vector<tuple<int,int,int,real_t> > const & contacts,
		vector<tuple<int, int, int> > const & obs,
		vector<tuple<int, Pi, Pr> > const & individuals) : params(params)
{
	Tinf = -1;
	for (auto it = contacts.begin(); it != contacts.end(); ++it) {
		int i,j,t;
		real_t lambda;
		tie(i,j,t,lambda) = *it;
		Tinf = max(Tinf, t + 1);
		add_contact(i, j, t, lambda);
	}

	vector<vector<int> > tobs(nodes.size());
	vector<vector<int> > sobs(nodes.size());
	for (auto it = obs.begin(); it != obs.end(); ++it) {
		int i,s,t;
		tie(i,s,t) = *it;
		Tinf = max(Tinf, t + 1);
		i = add_node(i);
		if (nodes.size() > tobs.size()) {
			tobs.resize(nodes.size());
			sobs.resize(nodes.size());
		}

		if (tobs[i].empty() || t >= tobs[i].back()) {
			tobs[i].push_back(t);
			sobs[i].push_back(s);
		} else {
			throw invalid_argument("time of observations should be ordered");
		}
	}

	for (auto it = individuals.begin(); it != individuals.end(); ++it) {
		int a = add_node(get<0>(*it));
		nodes[a].prob_i = get<1>(*it);
		nodes[a].prob_g = get<2>(*it);
	}

	for (int i = 0; i < int(nodes.size()); ++i) {
		vector<int> F = tobs[i];
		for (int k = 0; k < int(nodes[i].neighs.size()); ++k) {
			vector<int> const & tij = nodes[i].neighs[k].times;
			F.insert(F.end(), tij.begin(), tij.end());
		}
		sort(F.begin(), F.end());
		F.push_back(Tinf);
		nodes[i].times.push_back(-1);
		for (int k = 0; k < int(F.size()); ++k) {
			if (nodes[i].times.back() != F[k])
				nodes[i].times.push_back(F[k]);
		}
		int ntimes = nodes[i].times.size();
		nodes[i].bt.resize(ntimes, 1./ntimes);
		nodes[i].bg.resize(ntimes, 1./ntimes);
		nodes[i].ht.resize(ntimes);
		nodes[i].hg.resize(ntimes);
		set_field(i, tobs[i], sobs[i]);
		for (int k = 0; k  < int(nodes[i].neighs.size()); ++k) {
			omp_init_lock(&nodes[i].neighs[k].lock_);
			nodes[i].neighs[k].times.push_back(Tinf);
			nodes[i].neighs[k].lambdas.push_back(0.0);
			int nij = nodes[i].neighs[k].times.size();
			nodes[i].neighs[k].msg.resize(nij*nij);
		}
	}
	init();
	//showgraph();
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
	nodes.push_back(Node(i, params.prob_i, params.prob_r));
	return index[i];
}

void FactorGraph::add_contact(int i, int j, int t, real_t lambda)
{
	i = add_node(i);
	j = add_node(j);
	int ki = find_neighbor(i, j);
	int kj = find_neighbor(j, i);
	if (ki == int(nodes[i].neighs.size())) {
		assert(kj == int(nodes[j].neighs.size()));
		nodes[i].neighs.push_back(Neigh(j, kj));
		nodes[j].neighs.push_back(Neigh(i, ki));
	}
	Neigh & ni = nodes[i].neighs[ki];
	Neigh & nj = nodes[j].neighs[kj];
	if (ni.times.empty() || t > ni.times.back()) {
		ni.times.push_back(t);
		ni.lambdas.push_back(lambda);
		nj.times.push_back(t);
		nj.lambdas.push_back(0.0);
	} else if (t == ni.times.back()) {
		ni.lambdas.back() = lambda;
	} else {
		throw invalid_argument("time of contacts should be ordered");
	}
}

void FactorGraph::set_field(int i, vector<int> const & tobs, vector<int> const & sobs)
{
	// this assumes ordered observation times
	int tl = 0, gl = 0;
	int tu = nodes[i].times.size();
	int gu = nodes[i].times.size();
	int t = 0;
	for (int k = 0; k < int(tobs.size()); ++k) {
		int state = sobs[k];
		int to = tobs[k];
		while (nodes[i].times[t] != to)
			t++;
		switch (state) {
			case 0:
				tl = max(tl, t);
				gl = max(gl, t);
				break;
			case 1:
				tu = min(tu, t - 1);
				gl = max(gl, t);
				break;
			case 2:
				tu = min(tu, t - 1);
				gu = min(gu, t - 1);
				break;
			case -1:
				break;
		}
	}

	for(int t = 0; t < int(nodes[i].ht.size()); ++t) {
		nodes[i].ht[t] = (tl <= t && t <= tu);
		nodes[i].hg[t] = (gl <= t && t <= gu);
	}
}

void FactorGraph::show_graph()
{
	cerr << "Number of nodes " <<  int(nodes.size()) << endl;
	for(int i = 0; i < int(nodes.size()); i++) {
		cerr << "### index " << nodes[i].index << "###" << endl;
		cerr << "### in contact with " <<  int(nodes[i].neighs.size()) << "nodes" << endl;
		vector<Neigh> const & aux = nodes[i].neighs;
		for (int j = 0; j < int(aux.size()); j++) {
			cerr << "# neighbor " << nodes[aux[j].index].index << endl;
			cerr << "# in position " << aux[j].pos << endl;
			cerr << "# in contact " << int(aux[j].times.size()) << " times, in t: ";
			for (int t = 0; t < int(aux[j].times.size()); t++)
				cerr << aux[j].times[t] << " ";
			cerr << " " << endl;
		}
	}
}

void FactorGraph::show_beliefs(ostream & ofs)
{
	for(int i = 0; i < int(nodes.size()); ++i) {
		Node & f = nodes[i];
		ofs << "node " << f.index << ":" << endl;
		for (int t = 0; t < int(f.bt.size()); ++t) {
			ofs << "    " << f.times[t] << " " << f.bt[t] << " (" << f.ht[t] << ") " << f.bg[t] << " (" << f.hg[t] << ")" << endl;
		}
	}

}

void FactorGraph::show_msg(ostream & msgfile)
{
	for(int i = 0; i < int(nodes.size()); ++i) {
		for(int j = 0; j < int(nodes[i].neighs.size()); ++j) {
			for (int n = 0; n < int(nodes[i].neighs[j].msg.size()); ++n) {
				vector<real_t> & msg = nodes[i].neighs[j].msg;
				msgfile << msg[n] << " ";
			}
			msgfile << " " << endl;
		}
	}
}

void norm_msg(vector<real_t> & msg)
{
	real_t S = 0;
	for(int n = 0; n < int(msg.size()); ++n)
		S += msg[n];
	if (!(S > 0))
		throw domain_error("singularity error");
	for(int n = 0; n < int(msg.size()); ++n)
		msg[n] /= S;
}

real_t setmes(vector<real_t> & from, vector<real_t> & to, real_t damp)
{
	int n = from.size();
	real_t s = 0;
	for (int i = 0; i < n; ++i) {
		s += from[i];
	}
	real_t err = 0;
	for (int i = 0; i < n; ++i) {
        if (!(s > 0)){
            from[i] = 1./n;
            err = numeric_limits<real_t>::infinity();
        } else {
            from[i] /= s;
            err = max(err, abs(from[i] - to[i]));
        }
		to[i] = damp*to[i] + (1-damp)*from[i];
	}
	return err;
}

int Sij(Node const & f, Neigh const & v, int sij, int gi)
{
	// here gi stands for the ti + gi index
	return v.times[sij] <= f.times[gi] ? sij : v.times.size() - 1;
}

inline int idx(int sij, int sji, int qj)
{
	return sji + qj * sij;
}


ostream & operator<<(ostream & o, vector<real_t> const & m)
{
	o << "{";
	for (int i=0; i<int(m.size()); ++i){
		o << m[i] << " ";
	}
	o << "}";
	return o;
}

void FactorGraph::init()
{
	for(int i = 0; i < int(nodes.size()); ++i) {
		for(int j = 0; j < int(nodes[i].neighs.size()); ++j) {
			vector<real_t> & msg = nodes[i].neighs[j].msg;
			for(int ss = 0; ss < int(msg.size()); ++ss)
				msg[ss] = 1;
		}
		fill(nodes[i].bt.begin(), nodes[i].bt.end(), 1./nodes[i].bt.size());
		fill(nodes[i].bg.begin(), nodes[i].bg.end(), 1./nodes[i].bg.size());
	}
}

void update_limits(int ti, Node const &f, vector<int> & min_in, vector<int> & min_out)
{
	int n = min_in.size();
	for (int j = 0; j < n; ++j) {
		Neigh const & v = f.neighs[j];
		int const *b = &v.times[0];
		int qj = v.times.size();
		int const *e = &v.times[0] + qj;
		min_in[j] = min(qj - 1, int(std::lower_bound(b, e, f.times[ti]) - b));
		min_out[j] = min(qj - 1, int(std::upper_bound(b + min_in[j], e, f.times[ti]) - b));
	}
}

real_t FactorGraph::update(int i, real_t damping)
{
	Node & f = nodes[i];
	int const n = f.neighs.size();
	vector<vector<real_t> > UU(n);
	vector<vector<real_t> > HH(n);
	int const qi = f.bt.size();

	vector<real_t> ut(qi);
	vector<real_t> ug(qi);


	for (int j = 0; j < n; ++j) {
		Neigh & v = nodes[f.neighs[j].index].neighs[f.neighs[j].pos];
		omp_set_lock(&v.lock_);
		HH[j] = v.msg;
		omp_unset_lock(&v.lock_);
		UU[j].resize(v.msg.size());

	}
	vector<vector<real_t>> M = UU, R = UU;
	// proba tji >= ti for each j
	vector<real_t> C0(n), P0(n);
	// proba tji > ti for each j
	vector<real_t> C1(n), P1(n);

	vector<real_t> ht = f.ht;
	ht[0] *= params.pseed;
	for (int t = 1; t < qi - 1; ++t)
		ht[t] *= 1 - params.pseed - params.psus;
	ht[qi-1] *= params.psus;

	vector<int> min_in(n), min_out(n);
	real_t za = 0.0;

	for (int ti = 0; ti < qi; ++ti) if (f.ht[ti]) {
		update_limits(ti, f, min_in, min_out);

		for (int j = 0; j < n; ++j) {
			vector<real_t> & m = M[j];
			vector<real_t> & r = R[j];
			fill(m.begin(), m.end(), 0.0);
			fill(r.begin(), r.end(), 0.0);
			Neigh const & v = f.neighs[j];
			vector<real_t> const & h = HH[j];
			int const qj = v.times.size();
			for (int sji = min_in[j]; sji < qj; ++sji) {
				real_t pi = 1;
				for (int sij = min_out[j]; sij < qj - 1; ++sij) {
					real_t const l = v.lambdas[sij] * f.prob_i(v.times[sij]-f.times[ti]);
					m[idx(sji, sij, qj)] = l * pi * h[idx(sji, sij, qj)];
					r[idx(sji, sij, qj)] = l * pi * h[idx(sji, qj - 1, qj)];;
					pi *= 1 - l;
				}
				m[idx(sji, qj - 1, qj)] = pi * h[idx(sji, qj - 1, qj)];
				r[idx(sji, qj - 1, qj)] = pi * h[idx(sji, qj - 1, qj)];
			}
			// accumulate
			for (int sji = qj - 2; sji >=  min_in[j]; --sji) {
				for (int sij = qj - 1; sij >= min_out[j]; --sij) {
					m[idx(sji, sij, qj)] += m[idx(sji + 1, sij, qj)];
					r[idx(sji, sij, qj)] += r[idx(sji + 1, sij, qj)];
				}
			}

			for (int sji = qj - 1; sji >=  min_in[j]; --sji) {
				for (int sij = qj - 2; sij >= min_out[j]; --sij) {
					m[idx(sji, sij, qj)] += m[idx(sji, sij + 1, qj)];
					r[idx(sji, sij, qj)] += r[idx(sji, sij + 1, qj)];
				}
			}
		}

/*
             .-----min_out
             |   .-- min_out_g
   sij       v   v
   . . . . . . . . .
sji. . . . . . . . .
   . . . . . . . . .
   . . . . . a a b b <- min_in
   . . . . . a a b b
   . . . . . c c d d <- min_out
   . . . . . c c d d
   . . . . . c c d d
   . . . . . c c d d


C0 = a + c + b' + d'
C1 = c + d'
*/
		for (int gi = ti; gi < qi; ++gi) if (f.hg[gi]) {
			for (int j = 0; j < n; ++j) {
				vector<real_t> & m = M[j];
				vector<real_t> & r = R[j];
				Neigh const & v = f.neighs[j];
				int const qj = v.times.size();
				int const * b = &v.times[0];
				int const min_out_g = min(qj - 1, int(std::upper_bound(b + min_out[j], b + qj, f.times[gi]) - b));
				C0[j] = m[idx(min_in[j], min_out[j], qj)] - m[idx(min_in[j], min_out_g, qj)] + r[idx(min_in[j], min_out_g, qj)];
				C1[j] = m[idx(min_out[j], min_out[j], qj)] - m[idx(min_out[j], min_out_g, qj)] + r[idx(min_out[j], min_out_g, qj)];
			}

			real_t p0full = cavity(C0.begin(), C0.end(), P0.begin(), 1.0, multiplies<real_t>());
			real_t p1full = cavity(C1.begin(), C1.end(), P1.begin(), 1.0, multiplies<real_t>());

			//messages to ti, gi
			real_t const g_prob = f.prob_g(f.times[gi] - f.times[ti]) - (gi + 1 == qi ? 0.0 : f.prob_g(f.times[gi + 1] - f.times[ti]));
			real_t const a = g_prob  * (ti == 0 || ti == qi - 1 ? p0full : p0full - p1full);

			ug[gi] += ht[ti] * a;
			ut[ti] += f.hg[gi] * a;
			za += ht[ti] * f.hg[gi] * a;

			//messages to sij, sji
			for (int j = 0; j < n; ++j) {
				Neigh const & v = f.neighs[j];
				int const qj = v.times.size();
				real_t const p0 = P0[j];
				real_t const p01 = p0 - P1[j];
				for (int sji = min_in[j]; sji < qj; ++sji) {
					real_t pi = ht[ti] * f.hg[gi] * g_prob * (ti == 0 || v.times[sji] == f.times[ti] ? p0 : p01);
					for (int s = min_out[j]; s < qj - 1; ++s) {
						real_t & Uij = UU[j][idx(Sij(f, v, s, gi), sji, qj)];
						real_t const l = f.prob_i(v.times[s]-f.times[ti]) *  v.lambdas[s];
						Uij += pi * l;
						pi *= 1 - l;
					}
					UU[j][idx(Sij(f, v, qj - 1, gi), sji, qj)] += pi;
				}
			}
		}
	}
	f.f_ = -log(za);
	//apply external fields on t,h
	for (int t = 0; t < qi; ++t) {
		ut[t] *= ht[t];
		ug[t] *= f.hg[t];
	}
	//compute marginals on t,g
	real_t diff = max(setmes(ut, f.bt, damping), setmes(ug, f.bg, damping));
	for (int j = 0; j < n; ++j) {
		Neigh & v = f.neighs[j];
		omp_set_lock(&v.lock_);
		// diff = max(diff, setmes(UU[j], v.msg, damping));
		setmes(UU[j], v.msg, damping);
		omp_unset_lock(&v.lock_);

		real_t zj = 0; // z_{(sij,sji)}}
		int const qj = v.times.size();
		for (int sij = 0; sij < qj; ++sij) {
			for (int sji = 0; sji < qj; ++sji) {
				zj += HH[j][idx(sij, sji, qj)]*v.msg[idx(sji, sij, qj)];
			}
		}
		f.f_ += 0.5*log(zj); // half is cancelled by z_{a,(sij,sji)}
	}

	return diff;

}

real_t FactorGraph::iteration(real_t damping)
{
	int const N = nodes.size();
	real_t err = 0.0;
	vector<int> perm(N);
	for(int i = 0; i < N; ++i)
		perm[i] = i;
	random_shuffle(perm.begin(), perm.end());
#pragma omp parallel for reduction(max:err)
	for(int i = 0; i < N; ++i)
		err = max(err, update(perm[i], damping));
	return err;
}

real_t FactorGraph::loglikelihood() const
{
	real_t L = 0;
	for(auto nit = nodes.begin(), nend = nodes.end(); nit != nend; ++nit)
		L -= nit->f_;
	return L;
}

real_t FactorGraph::iterate(int maxit, real_t tol, real_t damping)
{
	real_t err = std::numeric_limits<real_t>::infinity();
	for (int it = 1; it <= maxit; ++it) {
		err = iteration(damping);
		cout << "it: " << it << " err: " << err << endl;
		if (err < tol)
			break;
	}
	return err;
}



ostream & operator<<(ostream & ost, FactorGraph const & f)
{
	int nasym = 0;
	int nedge = 0;
	int ncont = 0;
	for(auto nit = f.nodes.begin(); nit != f.nodes.end(); ++nit) {
		for (auto vit = nit->neighs.begin(); vit != nit->neighs.end(); ++vit) {
                        if (vit->index < nit->index)
                                continue;
			++nedge;
			ncont += vit->lambdas.size() - 1;
			if (vit->lambdas != f.nodes[vit->index].neighs[vit->pos].lambdas)
				++nasym;
		}
	}

	return ost << "FactorGraph\n"
                << "            nodes: " << f.nodes.size() << "\n"
		<< "            edges: " << nedge << " ("  << nasym <<  " asymmetric)\n"
		<< "    time contacts: " << ncont;
}

