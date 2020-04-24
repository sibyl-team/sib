// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <sstream>
#include "bp.h"
#include "params.h"

namespace py = pybind11;

using namespace std;


string showparams(Params const & p)
{
    return "sib.Params(mu=" + to_string(p.mu) +
        ",pseed=" + to_string(p.pseed) + ")";
}

string showfg(FactorGraph const & f) {
    return "sib.FactorGraph with " + to_string(f.nodes.size()) + " nodes";
}

PYBIND11_MODULE(_sib, m) {
    py::class_<FactorGraph>(m, "FactorGraph")
        .def(py::init<vector<tuple<int,int,int,real_t> >, vector<tuple<int,int,int> >, Params const &>(), py::arg("contacts"), py::arg("observations"), py::arg("params"))
        .def("update", &FactorGraph::iteration)
        .def("bt", &FactorGraph::get_tbeliefs)
        .def("bg", &FactorGraph::get_gbeliefs)
        .def("reset", &FactorGraph::init)
        .def("__repr__", &showfg)
        .def_readwrite("params", &FactorGraph::params);


    py::class_<Params>(m, "Params")
        .def(py::init<real_t, real_t>(), py::arg("mu") = 0.01, py::arg("pseed") = 0.01)
        .def_readwrite("mu", &Params::mu)
        .def_readwrite("pseed", &Params::pseed)
        .def("__repr__", &showparams);

}
