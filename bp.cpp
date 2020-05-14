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


int const Tinf = 1000000;
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



void FactorGraph::append_observation(int i, int s, int t)
{
	add_node(i);
	Node & n = nodes[i];
        if (t < n.times[n.times.size() - 2]) {
		cerr << t << " " << n.times[n.times.size() - 2] << endl;
                throw invalid_argument("observation time too small");
	} else if (t > n.times[n.times.size() - 2]) {
		n.push_back_time(t);
                // adjust infinite times
                for (int j = 0; j < int(n.neighs.size()); ++j) {
                        n.neighs[j].t.back() = n.times.size() - 1;
                }
        }
        int qi = n.times.size();
        int tobs = qi - 2;
	int tl = 0, gl = 0;
	int tu = qi;
	int gu = qi;
        switch (s) {
                case 0:
                        tl = max(tl, tobs);
                        gl = max(gl, tobs);
                        break;
                case 1:
                        tu = min(tu, tobs - 1);
                        gl = max(gl, tobs);
                        break;
                case 2:
                        tu = min(tu, tobs - 1);
                        gu = min(gu, tobs - 1);
                        break;
                case -1:
                        break;

        }

	for(int t = 0; t < qi; ++t) {
		n.ht[t] *= (tl <= t && t <= tu);
		n.hg[t] *= (gl <= t && t <= gu);
	}
}

Mes & operator++(Mes & msg)
{
	int qj = msg.qj;
	msg.qj++;
	msg.resize(msg.qj * msg.qj);


	for (int sij = qj - 1; sij >= 0; --sij) {
		for (int sji = qj - 1; sji >= 0; --sji) {
			msg(sji, sij) = msg[(qj + 1) * sij + sji];
		}
	}
	for (int s = 0; s < int(msg.qj); ++s) {
		msg(s, qj) = msg(s, qj - 1);
		msg(qj, s) = msg(qj - 1, s);
                msg(qj, qj) = msg(qj - 1, qj - 1);
	}
	return msg;
}

void FactorGraph::append_contact(int i, int j, int t, real_t lambdaij, real_t lambdaji)
{
        add_node(i);
        add_node(j);
	Node & fi = nodes[i];
	Node & fj = nodes[j];
	int qi = fi.times.size();
	int qj = fj.times.size();
	if (fi.times[qi - 2] > t || fj.times[qj - 2] > t)
		throw invalid_argument("time of contacts should be ordered");

	int ki = find_neighbor(i, j);
	int kj = find_neighbor(j, i);

	if (ki == int(fi.neighs.size())) {
		assert(kj == int(fj.neighs.size()));
		fi.neighs.push_back(Neigh(j, kj));
		fj.neighs.push_back(Neigh(i, ki));
	}

	Neigh & ni = fi.neighs[ki];
	Neigh & nj = fj.neighs[kj];
	if (fi.times[qi - 2] < t) {
		fi.push_back_time(t);
                ++qi;
	}
	if (fj.times[qj - 2] < t) {
		fj.push_back_time(t);
                ++qj;
	}
	if (ni.t.size() < 2 || ni.t[ni.t.size() - 2] < qi - 2) {
		ni.t.back() = qi - 2;
		nj.t.back() = qj - 2;
		ni.t.push_back(qi - 1);
		nj.t.push_back(qj - 1);
		ni.lambdas.back() = lambdaij;
		nj.lambdas.back() = lambdaji;
                ni.lambdas.push_back(0.0);
                nj.lambdas.push_back(0.0);
		++ni.msg;
		++nj.msg;
	} else if (ni.t[ni.t.size() - 2] == qi - 2) {
		ni.lambdas[ni.t.size() - 2] = lambdaij;
	} else {
		throw invalid_argument("time of contacts should be ordered");
	}
        // adjust infinite times
        for (int k = 0; k < int(fi.neighs.size()); ++k)
                fi.neighs[k].t.back() = qi - 1;
        for (int k = 0; k < int(fj.neighs.size()); ++k)
                fj.neighs[k].t.back() = qj - 1;
}


FactorGraph::FactorGraph(Params const & params,
		vector<tuple<int,int,int,real_t> > const & contacts2,
		vector<tuple<int, int, int> > const & obs2,
		vector<tuple<int, std::shared_ptr<Proba>, std::shared_ptr<Proba>>> const & individuals) :
	params(params)
{
	for (auto it = individuals.begin(); it != individuals.end(); ++it) {
		add_node(get<0>(*it));
		Node & n = nodes[get<0>(*it)];
		n.prob_i = get<1>(*it);
		n.prob_r = get<2>(*it);
	}
	vector<tuple<int,int,int,real_t> > contacts = contacts2;
	sort(contacts.begin(), contacts.end(), [](tuple<int,int,int,real_t> const & x,tuple<int,int,int,real_t> const & y) {return get<2>(x) < get<2>(y);});

	vector<tuple<int, int, int> > obs = obs2;
	sort(obs.begin(), obs.end(), [](tuple<int,int,int> const & x,tuple<int,int,int> const & y) {return get<2>(x) < get<2>(y);});

	auto ic = contacts.begin(), ec = contacts.end();
	auto io = obs.begin(), eo = obs.end();
	while (ic != ec || io != eo) {
		int tc = ic == ec ? Tinf : get<2>(*ic);
		int to = io == eo ? Tinf : get<2>(*io);
		if (tc < to) {
			// cerr << "appending contact" << get<0>(*ic) << " " <<  get<1>(*ic)<< " " <<  get<2>(*ic) << " " <<  get<3>(*ic) << endl;
			append_contact(get<0>(*ic), get<1>(*ic), get<2>(*ic), get<3>(*ic));
			ic++;
		} else {
			// cerr << "appending obs" << get<0>(*io) << " " <<  get<1>(*io)<< " " <<  get<2>(*io)  << endl;
			append_observation(get<0>(*io), get<1>(*io), get<2>(*io));
			io++;
		}
	}
}

int FactorGraph::find_neighbor(int i, int j) const
{
	int k = 0;
	for (; k < int(nodes[i].neighs.size()); ++k)
		if (j == nodes[i].neighs[k].index)
			break;
	return k;
}

void FactorGraph::add_node(int i)
{
	for (int j = nodes.size(); j < i + 1; ++j)
		nodes.push_back(Node(params.prob_i, params.prob_r));
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
		cerr << "### index " << i << "###" << endl;
		cerr << "### in contact with " <<  int(nodes[i].neighs.size()) << "nodes" << endl;
		vector<Neigh> const & aux = nodes[i].neighs;
		for (int j = 0; j < int(aux.size()); j++) {
			cerr << "# neighbor " << aux[j].index << endl;
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
		ofs << "node " << i << ":" << endl;
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
			fill(CG0[j].begin(), CG0[j].end(), 0.0);
			fill(CG01[j].begin(), CG01[j].end(), 0.0);
		}
		min_g = min_out;
		for (int gi = ti; gi < qi; ++gi) if (f.hg[gi]) {
			for (int j = 0; j < n; ++j) {
				Mes & m = M[j];
				Mes & r = R[j];
				Neigh const & v = f.neighs[j];
				int const qj = v.t.size();
				int const *b = &v.t[0];
				min_g[j] = upper_bound(b + min_g[j], b + qj - 1, gi) - b;

				/*
				              .-----min_out
				              |   .-- min_g
				      sij     v   v
				      . . . . . . . .
				   sji. . . . . . . .
				      . . . . . . . .
				      . . . . a a b b <- min_in
				      . . . . c c d d <- min_out
				      . . . . c c d d
				      . . . . c c d d
				      . . . . c c d d


				   C0 = a + c + b' + d' = (a + c + b + d) - (b + d) + (b' + d')
				   C1 = c + d'          = c + d           - d       + d'
				*/
				C0[j] = m(min_in[j], min_out[j]) - m(min_in[j], min_g[j]) + r(min_in[j], min_g[j]);
				C1[j] = m(min_out[j], min_out[j]) - m(min_out[j], min_g[j]) + r(min_out[j], min_g[j]);
			}

			real_t p0full = cavity(C0.begin(), C0.end(), P0.begin(), 1.0, multiplies<real_t>());
			real_t p1full = cavity(C1.begin(), C1.end(), P1.begin(), 1.0, multiplies<real_t>());

			//messages to ti, gi
			real_t const pg = prob_r(f.times[gi] - f.times[ti]) - (gi >= qi - 1 ? 0.0 : prob_r(f.times[gi + 1] - f.times[ti]));
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
	for(int i = 0; i < int(f.nodes.size()); ++i) {
		for (auto vit = f.nodes[i].neighs.begin(), vend = f.nodes[i].neighs.end(); vit != vend; ++vit) {
                        if (vit->index < i)
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

