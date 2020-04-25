// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <sstream>
#include <numeric>
#include <iterator>
#include "bp.h"

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


PYBIND11_MODULE(_sib, m) {
    py::class_<FactorGraph>(m, "FactorGraph")
        .def(py::init<vector<tuple<int,int,int,real_t> >,
                vector<tuple<int,int,int> >,
                Params const &>(),
                py::arg("contacts"),
                py::arg("observations"),
                py::arg("params"))
        .def("update", &FactorGraph::iteration)
        .def("bt", &FactorGraph::get_tbeliefs)
        .def("bg", &FactorGraph::get_gbeliefs)
        .def("marginals", &get_marginals)
        .def("reset", &FactorGraph::init)
        .def("times", &get_times)
        .def("__repr__", &show_fg)
        .def_readwrite("params", &FactorGraph::params);


    py::class_<Params>(m, "Params")
        .def(py::init<real_t, real_t>(),
                py::arg("mu") = 0.01,
                py::arg("pseed") = 0.01)
        .def_readwrite("mu", &Params::mu)
        .def_readwrite("pseed", &Params::pseed)
        .def("__repr__", &show_params);

}
