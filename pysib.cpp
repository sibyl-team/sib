// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
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

tuple<real_t, real_t, real_t> get_marginal_t(Node const & n, int t)
{
        return get_marginal(n)[t];
}


int get_index(FactorGraph & f, int i)
{
    auto it = f.index.find(i);
    if (it == f.index.end())
       throw py::key_error("key not found");
    return it->second;
}


PYBIND11_MODULE(_sib, m) {
    py::bind_vector<std::vector<real_t>>(m, "VectorReal");
    py::bind_vector<std::vector<int>>(m, "VectorInt");
    py::bind_vector<std::vector<Node>>(m, "VectorNode");
    //py::bind_vector<std::vector<tuple<real_t, real_t, real_t>>(m, "VectorTuple");


    py::class_<FactorGraph>(m, "FactorGraph")
        .def(py::init<Params const &,
                vector<tuple<int,int,int,real_t>>,
                vector<tuple<int,int,int>>,
                vector<tuple<int,Pi,Pr>> >(),
                py::arg("params"),
                py::arg("contacts"),
                py::arg("observations"),
                py::arg("individuals") = vector<tuple<int,Pi,Pr>>())
        .def("update", &FactorGraph::iteration)
        .def("loglikelihood", &FactorGraph::loglikelihood)
        .def("reset", &FactorGraph::init)
        .def("get_index", &get_index)
        .def("__repr__", &print<FactorGraph>)
        .def_readonly("nodes", &FactorGraph::nodes)
        .def_readonly("params", &FactorGraph::params);

    py::class_<Node>(m, "Node")
        .def("marginal", &get_marginal)
        .def("marginal_t", &get_marginal_t)
        .def_readwrite("ht", &Node::ht)
        .def_readwrite("hg", &Node::hg)
        .def_readonly("bt", &Node::bt)
        .def_readonly("bg", &Node::bg)
        .def_readonly("times", &Node::times)
        .def_readonly("index", &Node::index);

    py::class_<Uniform>(m, "Uniform")
        .def(py::init<real_t>(), py::arg("p") = 1.0)
        .def_readwrite("p", &Uniform::p)
        .def("__repr__", &print<Uniform>);

    py::class_<Exponential>(m, "Exponential")
        .def(py::init<real_t>(), py::arg("mu") = 0.1)
        .def_readwrite("mu", &Exponential::mu)
        .def("__repr__", &print<Exponential>);

    py::class_<Gamma>(m, "Gamma")
        .def(py::init<real_t, real_t>(), py::arg("k") = 1.0, py::arg("mu") = 0.1)
        .def_readwrite("k", &Gamma::k)
        .def_readwrite("mu", &Gamma::mu)
        .def("__repr__", &print<Gamma>);

    py::class_<Params>(m, "Params")
        .def(py::init<Pi, Pr, real_t, real_t>(),
                "Params class. prob_i and prob_r parameters are defaults.",
                py::arg("prob_i") = Pi(1.0),
                py::arg("prob_r") = Pr(1.0, 0.1),
                py::arg("pseed") = 0.01,
                py::arg("psus") = 0.5)
        .def_readwrite("prob_r", &Params::prob_r)
        .def_readwrite("prob_i", &Params::prob_i)
        .def_readwrite("pseed", &Params::pseed)
        .def_readwrite("psus", &Params::psus)
        .def("__repr__", &print<Params>);
    m.def("set_num_threads", &omp_set_num_threads);
}
