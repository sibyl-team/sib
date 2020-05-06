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

void cumsum(Mes & m, int a, int b)
{
	for (int sij = m.qj - 2; sij >= b; --sij)
		m(m.qj - 1, sij) += m(m.qj -1, sij + 1);
	for (int sji = m.qj - 2; sji >= a; --sji) {
		real_t r = m(sji, m.qj - 1);
		m(sji, m.qj - 1) += m(sji + 1, m.qj - 1);
		for (int sij = m.qj - 2; sij >= b; --sij) {
			r += m(sji, sij);
			m(sji, sij) = r + m(sji + 1, sij);
		}
	}
}


FactorGraph::FactorGraph(Params const & params,
		vector<tuple<int,int,int,real_t> > const & contacts,
		vector<tuple<int, int, int> > const & obs,
		vector<tuple<int, std::shared_ptr<Proba>, std::shared_ptr<Proba>> > const & individuals) : params(params)
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
		nodes[a].prob_r = get<2>(*it);
	}

	if (nodes.size() > tobs.size()) {
		tobs.resize(nodes.size());
		sobs.resize(nodes.size());
	}
	for (int i = 0; i < int(nodes.size()); ++i) {
		Node & f = nodes[i];
		int const n = f.neighs.size();
		vector<tuple<int,int>> F;
		for (auto it = tobs[i].begin(); it != tobs[i].end(); ++it)
			F.push_back(make_tuple(*it, -1));
		for (int j = 0; j < n; ++j) {
			Neigh & v = f.neighs[j];
			v.t.push_back(Tinf);
			v.lambdas.push_back(0);
			for (int s = 0; s < int(v.t.size()); ++s)
				F.push_back(make_tuple(v.t[s], j));
		}
		sort(F.begin(), F.end());
		F.push_back(make_tuple(Tinf,-1));
		f.times.push_back(-1);
		vector<int> pointer(n);
		for (auto it = F.begin(); it != F.end(); ++it) {
			int t, j;
			tie(t,j) = *it;
			if (f.times.back() != t)
				f.times.push_back(t);
			if (j != -1)
				f.neighs[j].t[pointer[j]++] = f.times.size() - 1;
		}
		int qi = f.times.size();
		f.bt.resize(qi, 1./qi);
		f.bg.resize(qi, 1./qi);
		f.ht.resize(qi);
		f.hg.resize(qi);
		set_field(i, tobs[i], sobs[i]);
		for (int j = 0; j  < n; ++j) {
			Neigh & v = f.neighs[j];
			omp_init_lock(&v.lock_);
			int nij = v.t.size();
			nodes[i].neighs[j].msg = Mes(nij);
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
	if (ni.t.empty() || t > ni.t.back()) {
		ni.t.push_back(t);
		ni.lambdas.push_back(lambda);
		nj.t.push_back(t);
		nj.lambdas.push_back(0.0);
	} else if (t == ni.t.back()) {
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
		nodes[i].ht[t] = params.softconstraint + (1-params.softconstraint)*(tl <= t && t <= tu);
		nodes[i].hg[t] = params.softconstraint + (1-params.softconstraint)*(gl <= t && t <= gu);
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
			cerr << "# in contact " << int(aux[j].t.size()) << " times, in t: ";
			for (int s = 0; s < int(aux[j].t.size()); s++)
				cerr << aux[j].t[s] << " ";
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
				msgfile << nodes[i].neighs[j].msg[n] << " ";
			}
			msgfile << " " << endl;
		}
	}
}

void norm_msg(Mes & msg)
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
			Mes & msg = nodes[i].neighs[j].msg;
			fill(msg.begin(), msg.end(), 1. / msg.size());
		}
		fill(nodes[i].bt.begin(), nodes[i].bt.end(), 1. / nodes[i].bt.size());
		fill(nodes[i].bg.begin(), nodes[i].bg.end(), 1. / nodes[i].bg.size());
	}
}

void update_limits(int ti, Node const &f, vector<int> & min_in, vector<int> & min_out)
{
	int n = min_in.size();
	for (int j = 0; j < n; ++j) {
		Neigh const & v = f.neighs[j];
		int qj = v.t.size();
		int const *b = &v.t[0];
		int const *e = &v.t[0] + qj - 1;
		min_in[j] = lower_bound(b + min_in[j], e, ti) - b;
		min_out[j] = min_in[j] + (v.t[min_in[j]] == ti && min_in[j] < qj - 1);
	}
}


real_t FactorGraph::update(int i, real_t damping)
{
	Node & f = nodes[i];
	Proba const & prob_i = *f.prob_i;
	Proba const & prob_r = *f.prob_r;
	int const n = f.neighs.size();
	int const qi = f.bt.size();

	// allocate buffers
	vector<Mes> UU, HH, M, R;
	vector<real_t> ut(qi);
	vector<real_t> ug(qi);
	for (int j = 0; j < n; ++j) {
		Neigh const & v = nodes[f.neighs[j].index].neighs[f.neighs[j].pos];
		v.lock();
		HH.push_back(v.msg);
		v.unlock();
		UU.push_back(Mes(v.t.size()));
		R.push_back(Mes(v.t.size()));
		M.push_back(Mes(v.t.size()));
	}
	vector<real_t> C0(n), P0(n); // probas tji >= ti for each j
	vector<real_t> C1(n), P1(n); // probas tji > ti for each j
	vector<vector<real_t>> CG0(n, vector<real_t>(qi));
	vector<vector<real_t>> CG01(n, vector<real_t>(qi));
	vector<int> min_in(n), min_out(n), min_g(n);
	vector<real_t> ht = f.ht;

	// apply external fields
	ht[0] *= params.pseed;
	for (int t = 1; t < qi - 1; ++t)
		ht[t] *= 1 - params.pseed - params.psus;
	ht[qi-1] *= params.psus;

	real_t za = 0.0;
	for (int ti = 0; ti < qi; ++ti) if (f.ht[ti]) {
		update_limits(ti, f, min_in, min_out);

		for (int j = 0; j < n; ++j) {
			Mes & m = M[j]; // no need to clear, just use the bottom right corner
			Mes & r = R[j];
			Neigh const & v = f.neighs[j];
			Mes const & h = HH[j];
			int const qj = h.qj;
			for (int sji = min_in[j]; sji < qj; ++sji) {
				real_t pi = 1;
				for (int sij = min_out[j]; sij < qj - 1; ++sij) {
					real_t const l =  prob_i(f.times[v.t[sij]]-f.times[ti], v.lambdas[sij]);
					m(sji, sij) = l * pi * h(sji, sij);
					r(sji, sij) = l * pi * h(sji, qj - 1);;
					pi *= 1 - l;
				}
				m(sji, qj - 1) = pi * h(sji, qj - 1);
				r(sji, qj - 1) = pi * h(sji, qj - 1);
			}
			cumsum(m, min_in[j], min_out[j]);
			cumsum(r, min_in[j], min_out[j]);
			fill(CG0[j].begin(), CG0[j].end(), real_t(0));
			fill(CG01[j].begin(), CG01[j].end(), real_t(0));
		}
		min_g = min_out;
		for (int gi = ti; gi < qi; ++gi) if (f.hg[gi]) {
			for (int j = 0; j < n; ++j) {
				Mes & m = M[j];
				Mes & r = R[j];
				Neigh const & v = f.neighs[j];
				int const qj = v.t.size();
				int const * b = &v.t[0];
				min_g[j] = upper_bound(b + min_g[j], b + qj - 1, gi) - b;

				/*
				                .-----min_out
				                |   .-- min_g
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


				   C0 = a + c + b' + d' = (a + c + b + d) - (b + d) + (b' + d')
				   C1 = c + d'          = c + d           - d       + d'
				*/
				C0[j] = m(min_in[j], min_out[j]) - m(min_in[j], min_g[j]) + r(min_in[j], min_g[j]);
				C1[j] = m(min_out[j], min_out[j]) - m(min_out[j], min_g[j]) + r(min_out[j], min_g[j]);
			}

			real_t p0full = cavity(C0.begin(), C0.end(), P0.begin(), 1.0, multiplies<real_t>());
			real_t p1full = cavity(C1.begin(), C1.end(), P1.begin(), 1.0, multiplies<real_t>());

			//messages to ti, gi
			real_t const pg = prob_r(f.times[gi] - f.times[ti]) - (gi + 1 == qi ? 0.0 : prob_r(f.times[gi + 1] - f.times[ti]));
			real_t const a = pg * (ti == 0 || ti == qi - 1 ? p0full : p0full - p1full);

			ug[gi] += ht[ti] * a;
			ut[ti] += f.hg[gi] * a;
			za += ht[ti] * f.hg[gi] * a;

			for (int j = 0; j < n; ++j) {
				CG0[j][gi] += P0[j] * ht[ti] * f.hg[gi] * pg;
				CG01[j][gi] += (P0[j]-P1[j]) * ht[ti] * f.hg[gi] * pg;
			}
		}
		//messages to sij, sji
		for (int j = 0; j < n; ++j) {
			partial_sum(CG0[j].rbegin(), CG0[j].rend(), CG0[j].rbegin());
			partial_sum(CG01[j].rbegin(), CG01[j].rend(), CG01[j].rbegin());
			Neigh const & v = f.neighs[j];
			int const qj = v.t.size();
			for (int sji = min_in[j]; sji < qj; ++sji) {
				// note: ti == qi - 1 implies ti == v.t[sji]
				vector<real_t> const & CG = ti == 0 || ti == v.t[sji] ? CG0[j] : CG01[j];
				real_t pi = 1;
				real_t c = 0;
				for (int sij = min_out[j]; sij < qj - 1; ++sij) {
					int const tij = v.t[sij];
					real_t const l = prob_i(f.times[tij] - f.times[ti],  v.lambdas[sij]);
					UU[j](sij, sji) += CG[tij] * pi * l;
					c += (CG[ti] - CG[tij]) * pi * l;
					pi *= 1 - l;
				}
				UU[j](qj - 1, sji) += c + CG[ti] * pi;
			}
		}
	}
	f.f_ = -log(za);
	//apply external fields on t,h
	for (int t = 0; t < qi; ++t) {
		ut[t] *= ht[t];
		ug[t] *= f.hg[t];
	}
	//compute beliefs on t,g
	real_t diff = max(setmes(ut, f.bt, damping), setmes(ug, f.bg, damping));
	for (int j = 0; j < n; ++j) {
		Neigh & v = f.neighs[j];
		v.lock();
		// diff = max(diff, setmes(UU[j], v.msg, damping));
		setmes(UU[j], v.msg, damping);
		v.unlock();

		real_t zj = 0; // z_{(sij,sji)}}
		int const qj = v.t.size();
		for (int sij = 0; sij < qj; ++sij) {
			for (int sji = 0; sji < qj; ++sji) {
				zj += HH[j](sij, sji)*v.msg(sji, sij);
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
	real_t err = numeric_limits<real_t>::infinity();
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

