// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <sstream>
#include <numeric>
#include <iterator>
#include <exception>
#include "bp.h"

#include <pybind11/stl_bind.h>

PYBIND11_MAKE_OPAQUE(std::vector<real_t>);


namespace py = pybind11;

using namespace std;


string show_params(Params const & p)
{
    return "sib.Params(mu=" + to_string(p.mu) +
        ",pseed=" + to_string(p.pseed) + ")";
}

string show_fg(FactorGraph const & f) {
    return "sib.FactorGraph with " + to_string(f.nodes.size()) + " nodes";
}

map<int, vector<int> >
get_times(FactorGraph const & f) {
    map<int, vector<int> > times;
    for (int i = 0; i < int(f.nodes.size()); ++i) {
        (times[f.nodes[i].index] = f.nodes[i].times).pop_back();
    }
    return times;
}

map<int, vector<tuple<real_t, real_t, real_t> > >
get_marginals(FactorGraph const & f)
{
    map<int, vector<tuple<real_t, real_t, real_t> > > marg;
    for (int i = 0; i < int(f.nodes.size()); ++i) {
        Node const & n = f.nodes[i];
        vector<real_t> rbt(n.bt.size());
        vector<real_t> lbg(n.bg.size());
        int const T = n.bt.size() - 1;
        lbg[0] = n.bg[0];
        rbt[T] = n.bt[T];
        for (int t = 1; t <= T; ++t) {
            lbg[t] = lbg[t-1] + n.bg[t];
            rbt[T - t] = rbt[T - t + 1] + n.bt[T - t];
        }
        marg[n.index] = vector<tuple<real_t, real_t, real_t>>(T + 1);
        marg[n.index][0] = make_tuple(rbt[0], 0.0, 1.0 - rbt[0]); // this is (1,0,0)
        for (int t = 1; t <= T; ++t)
            marg[n.index][t] = make_tuple(rbt[t], lbg[t-1], 1-rbt[t]-lbg[t-1]);
    }
    return marg;
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

    py::class_<FactorGraph>(m, "FactorGraph")
        .def(py::init<vector<tuple<int,int,int,real_t> >,
                vector<tuple<int,int,int> >,
                Params const &>(),
                py::arg("contacts"),
                py::arg("observations"),
                py::arg("params"))
        .def("update", &FactorGraph::iteration)
        .def("reset", &FactorGraph::init)
        .def("get_index", &get_index)
        .def("__repr__", &show_fg)
        .def_readwrite("nodes", &FactorGraph::nodes)
        .def_readonly("params", &FactorGraph::params);

    py::class_<Node>(m, "Node")
        .def_readwrite("ht", &Node::ht)
        .def_readwrite("hg", &Node::hg)
        .def_readonly("bt", &Node::bt)
        .def_readonly("bg", &Node::bg)
        .def_readonly("times", &Node::times)
        .def_readonly("index", &Node::index);

    py::class_<Params>(m, "Params")
        .def(py::init<real_t, real_t>(),
                py::arg("mu") = 0.01,
                py::arg("pseed") = 0.01)
        .def_readwrite("mu", &Params::mu)
        .def_readwrite("pseed", &Params::pseed)
        .def("__repr__", &show_params);

}
