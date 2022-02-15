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
#include <assert.h>
#include <tuple>
#include <exception>
#include "bp.h"
#include "cavity.h"

using namespace std;

int const Tinf = 1000000;

template <class T>
void cumsum(Message<T> &m, int a, int b)
{
	T r = m(0, 0);
	for (int sij = m.qj - 2; sij >= b; --sij)
		m(m.qj - 1, sij) += m(m.qj - 1, sij + 1);
	for (int sji = m.qj - 2; sji >= a; --sji)
	{
		r = m(sji, m.qj - 1);
		m(sji, m.qj - 1) += m(sji + 1, m.qj - 1);
		for (int sij = m.qj - 2; sij >= b; --sij)
		{
			r += m(sji, sij);
			m(sji, sij) = r + m(sji + 1, sij);
		}
	}
}

void FactorGraph::append_time(int i, times_t t)
{
	add_node(i);
	Node &n = nodes[i];
	// most common case
	if (t == n.times[n.times.size() - 2] || t == *lower_bound(n.times.begin(), n.times.end(), t))
		return;
	if (t > n.times[n.times.size() - 2])
	{
		n.push_back_time(t);
		// adjust infinite times
		for (int j = 0; j < int(n.neighs.size()); ++j)
		{
			n.neighs[j].t.back() = n.times.size() - 1;
		}
		return;
	}
	cerr << t << " < " << n.times[n.times.size() - 2] << endl;
	throw invalid_argument("observation time unexistent and too small");
}

void FactorGraph::append_observation(int i, shared_ptr<Test> const &o, times_t t)
{
	add_node(i);
	append_time(i, t);
	if (o != params.fakeobs)
		nodes[i].obs.push_back(make_tuple(t, o));
}

Mes &operator++(Mes &msg)
{
	int oldqj = msg.qj;
	msg.qj++;
	int qj = msg.qj;
	msg.resize(msg.qj * msg.qj);

	// msg(sji, sij) = msg[qj * sij + sji]
	for (int sij = oldqj - 1; sij >= 0; --sij)
	{
		for (int sji = oldqj - 1; sji >= 0; --sji)
		{
			msg(sji, sij) = msg[oldqj * sij + sji];
		}
	}
	msg(qj - 1, qj - 1) = msg(qj - 2, qj - 2);
	for (int s = 0; s < qj; ++s)
	{
		msg(s, qj - 1) = msg(s, qj - 2);
		msg(qj - 1, s) = msg(qj - 2, s);
	}
	return msg;
}

Mes &operator--(Mes &msg)
{
	int qj = msg.qj;
	msg.qj--;
	for (int sij = 0; sij < qj - 1; ++sij)
	{
		for (int sji = 0; sji < qj - 1; ++sji)
		{
			msg(sji, sij) = msg[qj * (sij + 1) + (sji + 1)];
		}
	}
	msg.resize(msg.qj * msg.qj);
	return msg;
}

void FactorGraph::drop_contacts(times_t t)
{
	for (size_t i = 0; i < nodes.size(); ++i)
	{
		Node &fi = nodes[i];
		for (size_t k = 0; k < fi.neighs.size(); ++k)
		{
			if (fi.times[fi.neighs[k].t[0]] < t)
				throw invalid_argument("can only drop first contact");
			else if (fi.times[fi.neighs[k].t[0]] == t)
			{
				fi.neighs[k].t.erase(fi.neighs[k].t.begin(), fi.neighs[k].t.begin() + 1);
				fi.neighs[k].lambdas.erase(fi.neighs[k].lambdas.begin(), fi.neighs[k].lambdas.begin() + 1);
				--fi.neighs[k].msg;
			}
		}
	}
}

void FactorGraph::append_contact(int i, int j, times_t t, real_t lambdaij, real_t lambdaji)
{
	if (i == j)
		throw invalid_argument("self loops are not allowed");
	add_node(i);
	add_node(j);
	Node &fi = nodes[i];
	Node &fj = nodes[j];
	int qi = fi.times.size();
	int qj = fj.times.size();
	if (fi.times[qi - 2] > t || fj.times[qj - 2] > t)
		throw invalid_argument("time of contacts should be ordered");

	int ki = find_neighbor(i, j);
	int kj = find_neighbor(j, i);

	if (ki == int(fi.neighs.size()))
	{
		assert(kj == int(fj.neighs.size()));
		fi.neighs.push_back(Neigh(j, kj));
		fj.neighs.push_back(Neigh(i, ki));
	}

	Neigh &ni = fi.neighs[ki];
	Neigh &nj = fj.neighs[kj];
	if (fi.times[qi - 2] < t)
	{
		fi.push_back_time(t);
		++qi;
	}
	if (fj.times[qj - 2] < t)
	{
		fj.push_back_time(t);
		++qj;
	}
	if (ni.t.size() < 2 || ni.t[ni.t.size() - 2] < qi - 2)
	{
		ni.t.back() = qi - 2;
		nj.t.back() = qj - 2;
		ni.t.push_back(qi - 1);
		nj.t.push_back(qj - 1);
		if (lambdaij != DO_NOT_OVERWRITE)
			ni.lambdas.back() = lambdaij;
		if (lambdaji != DO_NOT_OVERWRITE)
			nj.lambdas.back() = lambdaji;
		ni.lambdas.push_back(0.0);
		nj.lambdas.push_back(0.0);
		++ni.msg;
		++nj.msg;
	}
	else if (ni.t[ni.t.size() - 2] == qi - 2)
	{
		if (lambdaij != DO_NOT_OVERWRITE)
			ni.lambdas[ni.t.size() - 2] = lambdaij;
		if (lambdaji != DO_NOT_OVERWRITE)
			nj.lambdas[nj.t.size() - 2] = lambdaji;
	}
	else
	{
		throw invalid_argument("time of contacts should be ordered");
	}
	// adjust infinite times
	for (int k = 0; k < int(fi.neighs.size()); ++k)
	{
		fi.neighs[k].t.back() = qi - 1;
	}
	for (int k = 0; k < int(fj.neighs.size()); ++k)
	{
		fj.neighs[k].t.back() = qj - 1;
	}
}

FactorGraph::FactorGraph(Params const &params,
						 vector<tuple<int, int, times_t, real_t>> const &contacts,
						 vector<tuple<int, std::shared_ptr<Test>, times_t>> const &tests,
						 vector<tuple<int, std::shared_ptr<Proba>, std::shared_ptr<Proba>, std::shared_ptr<Proba>, std::shared_ptr<Proba>>> const &individuals) : params(params)
{
	for (auto it = individuals.begin(); it != individuals.end(); ++it)
	{
		if (!get<1>(*it) || !get<1>(*it) || !get<1>(*it) || !get<1>(*it))
			throw invalid_argument("invalid individual definition");
		add_node(get<0>(*it));
		Node &n = nodes[get<0>(*it)];
		n.prob_i = get<1>(*it);
		n.prob_r = get<2>(*it);
		n.prob_i0 = get<3>(*it);
		n.prob_r0 = get<4>(*it);
		n.df_i = RealParams(n.prob_i->theta.size());
		n.df_r = RealParams(n.prob_r->theta.size());
	}
	auto ic = contacts.begin(), ec = contacts.end();
	auto io = tests.begin(), eo = tests.end();
	while (ic != ec || io != eo)
	{
		int tc = ic == ec ? Tinf : get<2>(*ic);
		int to = io == eo ? Tinf : get<2>(*io);
		if (tc < to)
		{
			// cerr << "appending contact" << get<0>(*ic) << " " <<  get<1>(*ic)<< " " <<  get<2>(*ic) << " " <<  get<3>(*ic) << endl;
			append_contact(get<0>(*ic), get<1>(*ic), get<2>(*ic), get<3>(*ic));
			ic++;
		}
		else
		{
			// cerr << "appending obs" << get<0>(*io) << " " <<  get<1>(*io)<< " " <<  get<2>(*io)  << endl;
			append_time(get<0>(*io), get<2>(*io));
			io++;
		}
	}
	reset_observations(tests);
}

void FactorGraph::reset_observations(vector<tuple<int, shared_ptr<Test>, times_t>> const &obs)
{
	for (unsigned j = 0; j < nodes.size(); ++j)
		nodes[j].obs.clear();
	for (unsigned k = 0; k < obs.size(); ++k)
	{
		auto p = obs[k];
		int i = get<0>(p);
		auto o = get<1>(p);
		times_t t = get<2>(p);
		if (o != params.fakeobs)
			nodes[i].obs.push_back(make_tuple(t, o));
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
		nodes.push_back(Node(params.prob_i, params.prob_r, j));
}

void FactorGraph::show_graph()
{
	cerr << "Number of nodes " << int(nodes.size()) << endl;
	for (int i = 0; i < int(nodes.size()); i++)
	{
		cerr << "### index " << i << "###" << endl;
		cerr << "### in contact with " << int(nodes[i].neighs.size()) << "nodes" << endl;
		vector<Neigh> const &aux = nodes[i].neighs;
		for (int j = 0; j < int(aux.size()); j++)
		{
			cerr << "# neighbor " << aux[j].index << endl;
			cerr << "# in position " << aux[j].pos << endl;
			cerr << "# in contact " << int(aux[j].t.size()) << " times, in t: ";
			for (int s = 0; s < int(aux[j].t.size()); s++)
				cerr << aux[j].t[s] << " ";
			cerr << " " << endl;
		}
	}
}

void FactorGraph::show_beliefs(ostream &ofs)
{
	for (int i = 0; i < int(nodes.size()); ++i)
	{
		Node &f = nodes[i];
		ofs << "node " << i << ":" << endl;
		for (int t = 0; t < int(f.bt.size()); ++t)
		{
			ofs << "    " << f.times[t] << " " << f.bt[t] << " (" << f.ht[t] << ") " << f.bg[t] << " (" << f.hg[t] << ")" << endl;
		}
	}
}

void FactorGraph::show_msg(ostream &o)
{
	for (int i = 0; i < int(nodes.size()); ++i)
	{
		auto &n = nodes[i];
		for (int j = 0; j < int(n.neighs.size()); ++j)
		{
			auto &v = n.neighs[j];
			o << i << " <- " << v.index << " : " << endl;
			for (int sij = 0; sij < int(v.msg.qj); ++sij)
			{
				for (int sji = 0; sji < int(v.msg.qj); ++sji)
				{
					o << v.msg(sij, sji) << " ";
				}
				o << endl;
			}
		}
	}
}

void norm_msg(Mes &msg)
{
	real_t S = 0;
	for (int n = 0; n < int(msg.size()); ++n)
		S += msg[n];
	if (!(S > 0))
		throw domain_error("singularity error");
	for (int n = 0; n < int(msg.size()); ++n)
		msg[n] /= S;
}

real_t setmes(vector<real_t> &from, vector<real_t> &to, real_t damp)
{
	int n = from.size();
	real_t s = 0;
	for (int i = 0; i < n; ++i)
	{
		s += from[i];
	}
	real_t err = 0;
	for (int i = 0; i < n; ++i)
	{
		if (!(s > 0))
		{
			from[i] = 1. / n;
			err = numeric_limits<real_t>::infinity();
		}
		else
		{
			from[i] /= s;
			err = max(err, abs(from[i] - to[i]));
		}
		to[i] = damp * to[i] + (1 - damp) * from[i];
	}
	return err;
}

ostream &operator<<(ostream &o, vector<real_t> const &m)
{
	o << "{";
	for (int i = 0; i < int(m.size()); ++i)
	{
		o << m[i] << " ";
	}
	o << "}";
	return o;
}

void update_limits(int ti, Node const &f, vector<int> &min_in, vector<int> &min_out)
{
	int n = min_in.size();
	for (int j = 0; j < n; ++j)
	{
		Neigh const &v = f.neighs[j];
		int qj = v.t.size();
		int const *b = &v.t[0];
		int const *e = &v.t[0] + qj - 1;
		min_in[j] = lower_bound(b + min_in[j], e, ti) - b;
		min_out[j] = min_in[j] + (v.t[min_in[j]] == ti && min_in[j] < qj - 1);
	}
}

real_t FactorGraph::update(int i, real_t damping, bool learn)
{
	Node &f = nodes[i];
	auto const &obs = f.obs;
	int const n = f.neighs.size();
	int const qi = f.bt.size();

	RealParams const zero_r = RealParams(0.0, f.prob_r->theta.size());
	RealParams const zero_i = RealParams(0.0, f.prob_i->theta.size());
	// allocate buffers
	vector<Mes> UU, HH, M, R;
	vector<Message<RealParams>> dM, dR;
	vector<real_t> ut(qi), ug(qi);
	vector<vector<real_t>> CG0, CG01;
	vector<RealParams> dC0, dC1;

	for (int j = 0; j < n; ++j)
	{
		Neigh const &v = nodes[f.neighs[j].index].neighs[f.neighs[j].pos];
		v.lock();
		HH.push_back(v.msg);
		v.unlock();
		UU.push_back(Mes(v.t.size()));
		R.push_back(Mes(v.t.size()));
		M.push_back(Mes(v.t.size()));
		CG0.push_back(vector<real_t>(v.t.size() + 1));
		CG01.push_back(vector<real_t>(v.t.size() + 1));
		if (learn)
		{
			dR.push_back(Message<RealParams>(v.t.size(), zero_r));
			dM.push_back(Message<RealParams>(v.t.size(), zero_r));
			dC0.push_back(zero_i);
			dC1.push_back(zero_i);
		}
	}
	vector<real_t> C0(n), P0(n); // probas tji >= ti for each j
	vector<real_t> C1(n), P1(n); // probas tji > ti for each j
	vector<int> min_in(n), min_out(n);

	// main loop
	real_t za = 0.0;
	RealParams dzr = zero_r, dp1 = zero_r, dp2 = zero_r;
	RealParams dzi = zero_i, dl = zero_i, dpi = zero_i, dlpi = zero_i;
	real_t qauto = 1.0;
	for (int ti = 0; ti < qi; ++ti)
	{
		Proba const &prob_i = ti ? *f.prob_i : *f.prob_i0;
		Proba const &prob_r = ti ? *f.prob_r : *f.prob_r0;
		bool const dolearn = (ti > 0) && learn;
		real_t const pauto = (0 < ti && ti < qi - 1) ? params.pautoinf : 0.0;
		update_limits(ti, f, min_in, min_out);

		for (int j = 0; j < n; ++j)
		{
			Mes &m = M[j]; // no need to clear, just use the bottom right corner
			Mes &r = R[j];
			Neigh const &v = f.neighs[j];
			Mes const &h = HH[j];
			int const qj = h.qj;

			real_t pi = 1;
			dpi = zero_i;

			Message<RealParams> &dm = dM[j];
			Message<RealParams> &dr = dR[j];
			for (int sij = min_out[j]; sij < qj - 1; ++sij)
			{
				int tij = v.t[sij];
				real_t const l = prob_i(f.times[tij] - f.times[ti]) * v.lambdas[sij];
				for (int sji = min_in[j]; sji < qj; ++sji)
				{
					m(sji, sij) = l * pi * h(sji, sij);
					r(sji, sij) = l * pi * h(sji, qj - 1);
				}
				if (dolearn)
				{
					prob_i.grad(dl, f.times[tij] - f.times[ti]);
					dl *= v.lambdas[sij];
					dlpi = dl * pi + l * dpi;
					for (int sji = min_in[j]; sji < qj; ++sji)
					{
						// grad m & r
						dm(sji, sij) = dlpi * h(sji, sij);
						dr(sji, sij) = dlpi * h(sji, qj - 1);
					}
					dpi = dpi * (1 - l) - pi * dl;
				}
				pi *= 1 - l;
			}

			for (int sji = min_in[j]; sji < qj; ++sji)
			{
				m(sji, qj - 1) = pi * h(sji, qj - 1);
				r(sji, qj - 1) = pi * h(sji, qj - 1);
				if (dolearn)
				{
					dm(sji, qj - 1) = dpi * h(sji, qj - 1);
					dr(sji, qj - 1) = dpi * h(sji, qj - 1);
				}
			}

			cumsum(m, min_in[j], min_out[j]);
			cumsum(r, min_in[j], min_out[j]);
			// grad m & r
			if (dolearn)
			{
				cumsum(dm, min_in[j], min_out[j]);
				cumsum(dr, min_in[j], min_out[j]);
			}
			fill(CG01[j].begin(), CG01[j].end(), 0.0);
			fill(CG0[j].begin(), CG0[j].end(), 0.0);
		}
		auto min_g = min_out;
		real_t p0full = 0.0, p1full = 0.0;
		bool changed = true;
		for (int j = 0; j < n; ++j)
			--min_g[j];

		for (int gi = ti; gi < qi; ++gi)
		{
			real_t w = f.ht[ti] * f.hg[gi];
			if (ti == 0)
				w *= params.pseed;
			else if (ti == qi - 1)
				w *= params.psus;
			else
				w *= 1 - params.pseed - params.psus;

			for (unsigned k = 0; k < obs.size(); ++k)
			{
				times_t const t = get<0>(obs[k]);
				real_t const ps = get<1>(obs[k])->ps;
				real_t const pi = get<1>(obs[k])->pi;
				real_t const pr = get<1>(obs[k])->pr;
				w *= ps * (f.times[ti] >= t) + pi * (f.times[ti] < t && t <= f.times[gi]) + pr * (t > f.times[gi]);
			}
			for (int j = 0; j < n; ++j)
			{
				Neigh const &v = f.neighs[j];
				int const qj = v.t.size();
				int const *b = &v.t[0];
				int newming = upper_bound(b + max(0, min_g[j]), b + qj - 1, gi) - b;
				if (newming == min_g[j])
					continue;
				min_g[j] = newming;
				changed = true;
				Mes &m = M[j];
				Mes &r = R[j];
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
				// grad C
				if (dolearn)
				{
					auto &dm = dM[j];
					auto &dr = dR[j];
					dC0[j] = dm(min_in[j], min_out[j]) - dm(min_in[j], min_g[j]) + dr(min_in[j], min_g[j]);
					dC1[j] = dm(min_out[j], min_out[j]) - dm(min_out[j], min_g[j]) + dr(min_out[j], min_g[j]);
				}
			}
			if (changed)
			{
				changed = false;
				p0full = cavity(C0.begin(), C0.end(), P0.begin(), 1.0, multiplies<real_t>());
				p1full = cavity(C1.begin(), C1.end(), P1.begin(), 1.0, multiplies<real_t>());
			}
			// messages to ti, gi
			auto const d1 = f.times[gi] - f.times[ti];
			real_t const pg = gi < qi - 1 ? prob_r(d1) - prob_r(f.times[gi + 1] - f.times[ti]) : prob_r(d1);
			real_t const c = qauto * (ti == 0 || ti == qi - 1 ? p0full : (p0full - p1full * (1 - pauto)));
			real_t const b = w * pg;
			ug[gi] += b * c;
			ut[ti] += b * c;
			za += b * c;
			for (int j = 0; j < n; ++j)
			{
				CG0[j][min_g[j]] += b * qauto * P0[j];
				CG01[j][min_g[j]] += b * qauto * (P0[j] - P1[j] * (1 - pauto));
			}
			if (dolearn)
			{
				// grad theta_r
				prob_r.grad(dp1, d1);
				if (gi < qi - 1)
				{
					auto const d2 = f.times[gi + 1] - f.times[ti];
					prob_r.grad(dp2, d2);
					dzr += w * (dp1 - dp2) * c;
				}
				else
				{
					dzr += w * dp1 * c;
				}
				// grad theta_i
				for (int j = 0; j < n; ++j)
				{
					dzi += b * qauto * P0[j] * dC0[j];
					if (0 < ti && ti < qi - 1)
						dzi -= b * qauto * P1[j] * dC1[j] * (1 - pauto);
				}
			}
		}
		// messages to sij, sji
		for (int j = 0; j < n; ++j)
		{
			partial_sum(CG0[j].rbegin(), CG0[j].rend(), CG0[j].rbegin());
			partial_sum(CG01[j].rbegin(), CG01[j].rend(), CG01[j].rbegin());
			Neigh const &v = f.neighs[j];
			int const qj = v.t.size();
			for (int sji = min_in[j]; sji < qj; ++sji)
			{
				// note: ti == qi - 1 implies ti == v.t[sji]
				vector<real_t> const &CG = ti == 0 || ti == v.t[sji] ? CG0[j] : CG01[j];
				real_t pi = 1;
				real_t c = 0;
				for (int sij = min_out[j]; sij < qj - 1; ++sij)
				{
					int const tij = v.t[sij];
					real_t const l = prob_i(f.times[tij] - f.times[ti]) * v.lambdas[sij];
					// note: CG[sij + 1] counts everything with gi >= sij
					UU[j](sij, sji) += CG[sij + 1] * pi * l;
					c += (CG[0] - CG[sij + 1]) * pi * l;
					pi *= 1 - l;
				}
				UU[j](qj - 1, sji) += c + CG[0] * pi;
			}
		}
		qauto *= 1 - pauto;
	}
	f.f_ = log(za);
	// update parameters
	if (learn && za)
	{
		f.df_r = dzr / za;
		f.df_i = dzi / za;
	}

	// compute beliefs on t,g
	real_t diff = max(setmes(ut, f.bt, damping), setmes(ug, f.bg, damping));
	f.err_ = diff;
	for (int j = 0; j < n; ++j)
	{
		Neigh &v = f.neighs[j];
		v.lock();
		// diff = max(diff, setmes(UU[j], v.msg, damping));
		setmes(UU[j], v.msg, damping);
		v.unlock();

		real_t zj = 0; // z_{(sij,sji)}}
		int const qj = v.t.size();
		for (int sij = 0; sij < qj; ++sij)
		{
			for (int sji = 0; sji < qj; ++sji)
			{
				zj += HH[j](sij, sji) * v.msg(sji, sij);
			}
		}
		f.f_ -= 0.5 * log(zj); // half is cancelled by z_{a,(sij,sji)}
	}

	return diff;
}

real_t FactorGraph::iteration(real_t damping, bool learn)
{
	int const N = nodes.size();
	real_t err = 0.0;
	vector<int> perm(N);
	for (int i = 0; i < N; ++i)
		perm[i] = i;
	random_shuffle(perm.begin(), perm.end());
#pragma omp parallel for reduction(max \
								   : err)
	for (int i = 0; i < N; ++i)
		err = max(err, update(perm[i], damping, learn));
	return err;
}

real_t FactorGraph::loglikelihood() const
{
	real_t L = 0;
	for (auto nit = nodes.begin(), nend = nodes.end(); nit != nend; ++nit)
		L += nit->f_;
	return L;
}

real_t FactorGraph::iterate(int maxit, real_t tol, real_t damping, bool learn)
{
	real_t err = numeric_limits<real_t>::infinity();
	for (int it = 1; it <= maxit; ++it)
	{
		err = iteration(damping, learn);
		cout << "it: " << it << " err: " << err << endl;
		if (err < tol)
			break;
	}
	return err;
}

ostream &operator<<(ostream &ost, FactorGraph const &f)
{
	int nasym = 0;
	int nedge = 0;
	int ncont = 0;
	for (int i = 0; i < int(f.nodes.size()); ++i)
	{
		for (auto vit = f.nodes[i].neighs.begin(), vend = f.nodes[i].neighs.end(); vit != vend; ++vit)
		{
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
			   << "            edges: " << nedge << " (" << nasym << " asymmetric)\n"
			   << "    time contacts: " << ncont;
}
