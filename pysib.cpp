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
        ",pseed=" + to_string(p.pseed) +
        ",tol=" + to_string(p.tol) +
        ",maxit=" + to_string(p.maxit) + ")";
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
        .def(py::init<real_t, real_t, real_t, int>(), py::arg("mu") = 0.01, py::arg("pseed") = 0.01, py::arg("tol") = 1e-5, py::arg("maxit")=100)
        .def_readwrite("mu", &Params::mu)
        .def_readwrite("tol", &Params::tol)
        .def_readwrite("pseed", &Params::pseed)
        .def_readwrite("maxit", &Params::maxit)
        .def("__repr__", &showparams);

}
