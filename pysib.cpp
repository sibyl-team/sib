#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "bp.h"
#include "params.h"

namespace py = pybind11;

using namespace std;

PYBIND11_MODULE(sib, m) {
    py::class_<FactorGraph>(m, "FactorGraph")
        .def(py::init<vector<tuple<int,int,int,real_t> >, vector<tuple<int,int,int> >, Params const &>(), py::arg("contacts"), py::arg("observations"), py::arg("params"))
        .def("iterate", &FactorGraph::iterate)
        .def("bt", &FactorGraph::get_tbeliefs)
        .def("bg", &FactorGraph::get_gbeliefs);


    py::class_<Params>(m, "Params")
        .def(py::init<real_t, real_t, real_t, int>(), py::arg("mu"), py::arg("pseed"), py::arg("tol"), py::arg("maxit"));
}
