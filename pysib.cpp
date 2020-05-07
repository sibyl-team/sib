// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <string>
#include <sstream>
#include <numeric>
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <exception>
#include "bp.h"


PYBIND11_MAKE_OPAQUE(std::vector<real_t>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<Node>);
//PYBIND11_MAKE_OPAQUE(std::vector<tuple<real_t, real_t, real_t>>);

namespace py = pybind11;
using namespace std;
using boost::lexical_cast;

vector<real_t> make_vector(py::list & l)
{
    vector<real_t> v(l.size());
    int i = 0;
    for (py::handle o : l) {
        v[i++] = py::cast<real_t>(o);
    }
    return v;
}

PriorDiscrete make_discrete(py::list & l)
{
    return PriorDiscrete(make_vector(l));
}


ExpDiscrete make_exp_discrete(py::list & l)
{
    return ExpDiscrete(make_vector(l));
}



template<class T> string print(const T & t) { return lexical_cast<string>(t); }

map<int, vector<int> >
get_times(FactorGraph const & f) {
    map<int, vector<int> > times;
    for (int i = 0; i < int(f.nodes.size()); ++i) {
        (times[f.nodes[i].index] = f.nodes[i].times).pop_back();
    }
    return times;
}

vector<tuple<real_t, real_t, real_t>>
get_marginal(Node const & n)
{
        vector<real_t> rbt(n.bt.size());
        vector<real_t> lbg(n.bg.size());
        int const T = n.bt.size() - 1;
        lbg[0] = n.bg[0];
        rbt[T] = n.bt[T];
        for (int t = 1; t <= T; ++t) {
            lbg[t] = lbg[t-1] + n.bg[t];
            rbt[T - t] = rbt[T - t + 1] + n.bt[T - t];
        }
        auto marg = vector<tuple<real_t, real_t, real_t>>(T-1);
        for (int t = 1; t < T; ++t)
            marg[t-1] = make_tuple(rbt[t], 1-rbt[t]-lbg[t-1], lbg[t-1]);
        return marg;
}

tuple<real_t, real_t, real_t> get_marginal_index(Node const & n, int t)
{
        if (t < 0 || t >= int(n.bt.size()) - 2)
            throw py::key_error("Time out of range");
        else
            return get_marginal(n)[t];
}


int get_index(FactorGraph & f, int i)
{
    auto it = f.index.find(i);
    if (it == f.index.end())
       throw py::key_error("key not found");
    return it->second;
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
		msg(s, qj) = msg(s, qj -1);
		msg(s, qj - 1) = 0.0;
		msg(qj, s) = msg(qj - 1, s);
		msg(s, qj - 1) = 0.0;
	}
	return msg;
}

void append_observation(FactorGraph & G, int i, int s, int t)
{
        auto mi = G.index.find(i);
        if (mi == G.index.end())
                throw invalid_argument("unexistant node");

        i = mi->second;
        Node & n = G.nodes[i];
        n.times.back() = t;
        n.times.push_back(G.Tinf);
        t = n.bt.size();
        n.ht.push_back(1.0);
        n.hg.push_back(1.0);
        n.bt.push_back(1.0);
        n.bg.push_back(1.0);
        int qi = n.times.size();
	int tl = 0, gl = 0;
	int tu = qi;
	int gu = qi;
        switch(s) {
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
        fill(&n.ht[0], &n.ht[0] + tl, 0.0);
        fill(&n.ht[0] + tu, &n.ht[0] + qi, 0.0);
        fill(&n.bt[0], &n.bt[0] + gl, 0.0);
        fill(&n.bt[0] + gu, &n.bt[0] + qi, 0.0);
}

void append_contact(FactorGraph & G, int i, int j, int t, real_t lambda)
{
        G.Tinf = max(G.Tinf, t + 1);
	auto mi = G.index.find(i);
	auto mj = G.index.find(j);
	if (mi == G.index.end() || mj == G.index.end())
		throw invalid_argument("cannot append contact to unexistant node");
	i = mi->second;
	j = mj->second;
	Node & fi = G.nodes[i];
	Node & fj = G.nodes[j];
	int qi = fi.times.size();
	int qj = fj.times.size();
	if (fi.times[qi - 2] > t || fj.times[qj - 2] > t)
		throw invalid_argument("time of contacts should be ordered");
        G.Tinf = max(G.Tinf, t + 1);

	int ki = G.find_neighbor(i, j);
	int kj = G.find_neighbor(j, i);

	if (ki == int(fi.neighs.size())) {
		assert(kj == int(fj.neighs.size()));
		fi.neighs.push_back(Neigh(j, kj));
		fj.neighs.push_back(Neigh(i, ki));
                fi.neighs.back().msg = Mes(1);
                fj.neighs.back().msg = Mes(1);
                fi.neighs.back().msg[0] = 1.0;
                fj.neighs.back().msg[0] = 1.0;
		fi.neighs.back().t.push_back(G.Tinf);
		fj.neighs.back().t.push_back(G.Tinf);
                fi.neighs.back().lambdas.push_back(0);
                fj.neighs.back().lambdas.push_back(0);
	}

	Neigh & ni = fi.neighs[ki];
	Neigh & nj = fj.neighs[kj];
	if (fi.times[qi - 2] < t) {
		fi.times.back() = t;
		fi.times.push_back(G.Tinf);
                fi.ht.push_back(0);
                fi.hg.push_back(0);
                fi.bt.push_back(0);
                fi.bg.push_back(0);
                ++qi;
	}
	if (fj.times[qj - 2] < t) {
		fj.times.back() = t;
		fj.times.push_back(G.Tinf);
                fj.ht.push_back(0);
                fj.hg.push_back(0);
                fj.bt.push_back(0);
                fj.bg.push_back(0);
                ++qj;
	}
	if (ni.t[ni.t.size() - 2] < qi - 2) {
		ni.t.back() = qi - 2;
		nj.t.back() = qj - 2;
		ni.t.push_back(qi - 1);
		nj.t.push_back(qj - 1);
		ni.lambdas.back() = lambda;
		nj.lambdas.back() = 0.0;
                ni.lambdas.push_back(0.0);
                nj.lambdas.push_back(0.0);
		++ni.msg;
		++nj.msg;
	} else if (ni.t[ni.t.size() - 2] == qi - 2) {
		ni.lambdas[ni.t.size() - 2] = lambda;
	} else {
		throw invalid_argument("time of contacts should be ordered");
	}
}







PYBIND11_MODULE(_sib, m) {
    py::bind_vector<std::vector<real_t>>(m, "VectorReal");
    py::bind_vector<std::vector<int>>(m, "VectorInt");
    py::bind_vector<std::vector<Node>>(m, "VectorNode");
    //py::bind_vector<std::vector<tuple<real_t, real_t, real_t>>(m, "VectorTuple");

    py::class_<Proba, shared_ptr<Proba>>(m, "Proba");

    py::class_<Uniform, Proba, shared_ptr<Uniform>>(m, "Uniform")
        .def(py::init<real_t>(), py::arg("p") = 1.0)
        .def_readwrite("p", &Uniform::p)
        .def("__repr__", &print<Uniform>);

    py::class_<Exponential, Proba, shared_ptr<Exponential>>(m, "Exponential")
        .def(py::init<real_t>(), py::arg("mu") = 0.1)
        .def_readwrite("mu", &Exponential::mu)
        .def("__repr__", &print<Exponential>);

    py::class_<Gamma, Proba, shared_ptr<Gamma>>(m, "Gamma")
        .def(py::init<real_t, real_t>(), py::arg("k") = 1.0, py::arg("mu") = 0.1)
        .def_readwrite("k", &Gamma::k)
        .def_readwrite("mu", &Gamma::mu)
        .def("__repr__", &print<Gamma>);

    py::class_<PriorDiscrete, Proba, shared_ptr<PriorDiscrete>>(m, "PriorDiscrete")
        .def(py::init(&make_discrete))
        .def_readwrite("p", &PriorDiscrete::p)
        .def("__repr__", &print<PriorDiscrete>);

    py::class_<ExpDiscrete, Proba, shared_ptr<ExpDiscrete>>(m, "ExpDiscrete")
        .def(py::init(&make_exp_discrete))
        .def_readwrite("p", &ExpDiscrete::p)
        .def("__repr__", &print<ExpDiscrete>);

    py::class_<Params>(m, "Params")
        .def(py::init<shared_ptr<Proba> const &, shared_ptr<Proba> const &, real_t, real_t, real_t>(),
                "Params class. prob_i and prob_r parameters are defaults.",
                py::arg("prob_i") = *new Uniform(1.0),
                py::arg("prob_r") = *new Exponential(0.1),
                py::arg("pseed") = 0.01,
                py::arg("psus") = 0.5,
                py::arg("softconstraint") = 0.0)
        .def_readwrite("prob_r", &Params::prob_r)
        .def_readwrite("prob_i", &Params::prob_i)
        .def_readwrite("pseed", &Params::pseed)
        .def_readwrite("psus", &Params::psus)
        .def_readwrite("softconstraint", &Params::softconstraint)
        .def("__repr__", &print<Params>);

    py::class_<FactorGraph>(m, "FactorGraph")
        .def(py::init<Params const &,
                vector<tuple<int,int,int,real_t>>,
                vector<tuple<int,int,int>>,
                vector<tuple<int,shared_ptr<Proba>,shared_ptr<Proba>>>
                >(),
                py::arg("params") = Params(shared_ptr<Proba>(new Uniform(1.0)), shared_ptr<Proba>(new Exponential(0.5)), 0.1, 0.45, 0.0),
                py::arg("contacts") = vector<tuple<int,int,int,real_t>>(),
                py::arg("observations") = vector<tuple<int,int,int>>(),
                py::arg("individuals") = vector<tuple<int,shared_ptr<Proba>,shared_ptr<Proba>>>())
        .def("update", &FactorGraph::iteration)
        .def("loglikelihood", &FactorGraph::loglikelihood)
        .def("reset", &FactorGraph::init)
        .def("get_index", &get_index)
        .def("__repr__", &print<FactorGraph>)
        .def("append_contact", &append_contact)
        .def("append_observation", &append_observation)
        .def_readonly("nodes", &FactorGraph::nodes)
        .def_readonly("params", &FactorGraph::params);
    py::class_<Node>(m, "Node")
        .def("marginal", &get_marginal)
        .def("marginal_index", &get_marginal_index)
        .def_readwrite("ht", &Node::ht)
        .def_readwrite("hg", &Node::hg)
        .def_readonly("bt", &Node::bt)
        .def_readonly("bg", &Node::bg)
        .def_readonly("times", &Node::times)
        .def_readonly("prob_i", &Node::prob_i)
        .def_readonly("prob_r", &Node::prob_r)
        .def_readonly("index", &Node::index);

    m.def("set_num_threads", &omp_set_num_threads);
}
