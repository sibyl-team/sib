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
            times[i] = f.nodes[i].times;
            times[i].pop_back();
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
        return get_marginal(n)[t];
}



void check_index(FactorGraph const & G, int i)
{
        if (i < 0 || i >= int(G.nodes.size()))
                throw invalid_argument("unexistent index");

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
                "SIB Params class. prob_i and prob_r parameters are defaults.",
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

    py::class_<FactorGraph>(m, "FactorGraph", "SIB class representing the graphical model of the epidemics")
        .def(py::init<Params const &,
                vector<tuple<int,int,int,real_t>>,
                vector<tuple<int,int,int>>,
                vector<tuple<int,shared_ptr<Proba>,shared_ptr<Proba>>>
                >(),
                py::arg("params") = Params(shared_ptr<Proba>(new Uniform(1.0)), shared_ptr<Proba>(new Exponential(0.5)), 0.1, 0.45, 0.0),
                py::arg("contacts") = vector<tuple<int,int,int,real_t>>(),
                py::arg("observations") = vector<tuple<int,int,int>>(),
                py::arg("individuals") = vector<tuple<int,shared_ptr<Proba>,shared_ptr<Proba>>>())
        .def("update", &FactorGraph::iteration, "perform one iteration")
        .def("loglikelihood", &FactorGraph::loglikelihood, "compute the bethe log-likelihood")
        .def("__repr__", &print<FactorGraph>)
        .def("append_contact", (void (FactorGraph::*)(int,int,int,real_t,real_t)) &FactorGraph::append_contact, "append contact (i,j,t,lambdaij,lambdaji)")
        .def("append_contact", (void (FactorGraph::*)(int,int,int,real_t)) &FactorGraph::append_contact, "append contact (i,j,t,lambdaij)")
        .def("append_observation", &FactorGraph::append_observation, "append an observation (i,state,t)")
        .def_readonly("nodes", &FactorGraph::nodes, "all nodes in this FactorGraph")
        .def_readonly("params", &FactorGraph::params, "parameters");
    py::class_<Node>(m, "Node", "SIB class representing an individual")
        .def("marginal", &get_marginal, "compute marginal probabilities (pS,pI,pR) corresponding to times n.times[1:]")
        .def("marginal_index", &get_marginal_index, "marginal at a given time (excluding time -1)")
        .def_readwrite("ht", &Node::ht, "external prior on ti")
        .def_readwrite("hg", &Node::hg, "external prior on gi")
        .def_readonly("bt", &Node::bt, "belief on ti")
        .def_readonly("bg", &Node::bg, "belief on gi")
        .def_readonly("times", &Node::times, "event times of this node")
        .def_readonly("index", &Node::index, "node index (deprecated, do not use)")
        .def_readonly("prob_i", &Node::prob_i, "probability of infection as function of t-ti")
        .def_readonly("prob_r", &Node::prob_r, "cumulative probability of recovery P(tr>t)");

    m.def("set_num_threads", &omp_set_num_threads, "sets the maximum number of simultaneous cpu threads");
}
