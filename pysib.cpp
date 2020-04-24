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
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<Node>);


namespace py = pybind11;

using namespace std;


string show_params(Params const & p)
{
    return "sib.Params(mu=" + to_string(p.mu) +
        ",pseed=" + to_string(p.pseed) +
        ",damping=" + to_string(p.damping) + ")";

}

string show_fg(FactorGraph const & f)
{
	int nasym = 0;
	int nedge = 0;
	int ncont = 0;
	for(auto nit = f.nodes.begin(); nit != f.nodes.end(); ++nit) {
		for (auto vit = nit->neighs.begin(); vit != nit->neighs.end(); ++vit) {
                        if (vit->index < nit->index)
                                continue;
			++nedge;
			ncont += vit->lambdas.size() - 1;
			if (vit->lambdas != f.nodes[vit->index].neighs[vit->pos].lambdas)
				++nasym;
		}
	}

	return "sib.FactorGraph\n"
                  "            nodes: " + to_string(f.nodes.size()) + "\n"
		+ "            edges: " + to_string(nedge) + " (" + to_string(nasym) + " assymetric)\n"
		+ "    time contacts: " + to_string(ncont);
}

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
            marg[t-1] = make_tuple(rbt[t], lbg[t-1], 1-rbt[t]-lbg[t-1]);
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
    py::bind_vector<std::vector<int>>(m, "VectorInt");
    py::bind_vector<std::vector<Node>>(m, "VectorNode");

    py::class_<FactorGraph>(m, "FactorGraph")
        .def(py::init<Params const &,
                vector<tuple<int,int,int,real_t>>,
                vector<tuple<int,int,int>>,
                vector<tuple<int,real_t,real_t>> >(),
                py::arg("params"),
                py::arg("contacts"),
                py::arg("observations"),
                py::arg("individuals") = vector<tuple<int,real_t,real_t>>())
        .def("update", &FactorGraph::iteration)
        .def("loglikelihood", &FactorGraph::loglikelihood)
        .def("reset", &FactorGraph::init)
        .def("get_index", &get_index)
        .def("__repr__", &show_fg)
        .def_readwrite("nodes", &FactorGraph::nodes)
        .def_readonly("params", &FactorGraph::params);

    py::class_<Node>(m, "Node")
        .def("marginal", &get_marginal)
        .def_readwrite("ht", &Node::ht)
        .def_readwrite("hg", &Node::hg)
        .def_readonly("bt", &Node::bt)
        .def_readonly("bg", &Node::bg)
        .def_readonly("times", &Node::times)
        .def_readonly("index", &Node::index);

    py::class_<Params>(m, "Params")
        .def(py::init<real_t, real_t, real_t, real_t>(),
                "Params class. mu and k parameters are defaults.",
                py::arg("mu") = 0.01,
                py::arg("k") = 1.0,
                py::arg("pseed") = 0.01,
                py::arg("damping") = 0.0)
        .def_readwrite("mu", &Params::mu)
        .def_readwrite("pseed", &Params::pseed)
        .def_readwrite("damping", &Params::damping)
        .def("__repr__", &show_params);
    m.def("set_num_threads", &omp_set_num_threads);
}
