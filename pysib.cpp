// This file is part of sibilla : inference in epidemics with Belief Propagation
// Author: Alfredo Braunstein


#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/pytypes.h>
#include <pybind11/iostream.h>
#include <string>
#include <sstream>
#include <numeric>
#include <boost/lexical_cast.hpp>
#include <iterator>
#include <exception>
#include "bp.h"
#include "drop.h"


PYBIND11_MAKE_OPAQUE(std::valarray<real_t>);
PYBIND11_MAKE_OPAQUE(std::vector<real_t>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<Node>);

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

map<int, vector<times_t> >
get_times(FactorGraph const & f) {
        map<int, vector<times_t> > times;
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


template<int i>
real_t mygetter(Proba & p)
{
    return p.theta[i];
}

template<int i>
void mysetter(Proba & p, real_t x)
{
    p.theta[i] = x;
}

PYBIND11_MODULE(_sib, m) {
    py::class_<RealParams>(m, "RealParams", py::buffer_protocol())
        .def(py::init([](py::buffer const b) {
                py::buffer_info info = b.request();
                if (info.format != py::format_descriptor<real_t>::format() || info.ndim != 1)
                throw std::runtime_error("Incompatible buffer format!");

                auto v = new RealParams(info.shape[0]);
                memcpy(&(*v)[0], info.ptr, sizeof(real_t) * (size_t) (v->size()));
                return v;
                }))
        .def(py::init([](vector<real_t> const & p)->RealParams {return RealParams(&p[0], p.size());}))
        .def(py::init([](py::list & l)->RealParams {auto v = make_vector(l); return RealParams(&v[0], v.size());}))
        .def_buffer([](RealParams &p) -> py::buffer_info {
            return py::buffer_info(
                &p[0],                               /* Pointer to buffer */
                sizeof(real_t),                          /* Size of one scalar */
                py::format_descriptor<real_t>::format(), /* Python struct-style format descriptor */
                1,                                      /* Number of dimensions */
                { p.size() },                 /* Buffer dimensions */
                { sizeof(real_t) }             /* Strides (in bytes) for each index */
                );
        })
        .def("__add__", [](RealParams & p, RealParams & q)->RealParams { return p + q; })
        .def("__getitem__", [](const RealParams &p, ssize_t i) {
                if (i > int(p.size()))
                    throw py::index_error();
                return p[i];
                })
        .def("__setitem__", [](RealParams &p, ssize_t i, real_t v) {
                if (i > int(p.size()))
                    throw py::index_error();
                p[i] = v;
                })
        .def("__repr__", [](RealParams &p) {
                    string s = "RealParams([";
                    for (size_t i = 0; i < p.size(); ++i)
                        s += (i ? ",":"") + lexical_cast<string>(p[i]);
                    s+="])";
                    return s;
                });
    // py::add_ostream_redirect(m, "ostream_redirect");
    py::bind_vector<std::vector<int>>(m, "VectorInt");
    py::bind_vector<std::vector<real_t>>(m, "VectorReal");
    py::bind_vector<std::vector<Node>>(m, "VectorNode");

    py::class_<Proba, shared_ptr<Proba>>(m, "Proba")
        .def("__call__", [](Proba const & p, real_t d) { return p(d); } )
        .def("grad", [](Proba const & p, real_t d) { RealParams dtheta(0.0, p.theta.size()); p.grad(dtheta, d); return dtheta;} )
        .def("__repr__", &lexical_cast<string, Proba>)
        .def_property("theta", &Proba::get_theta, &Proba::set_theta);

    py::class_<Uniform, Proba, shared_ptr<Uniform>>(m, "Uniform")
        .def(py::init<real_t>(), py::arg("p") = 1.0)
        .def_property("p", &mygetter<0>, &mysetter<0>);

    py::class_<Exponential, Proba, shared_ptr<Exponential>>(m, "Exponential")
        .def(py::init<real_t>(), py::arg("mu") = 0.1)
        .def_property("mu", &mygetter<0>, &mysetter<0>);

    py::class_<Gamma, Proba, shared_ptr<Gamma>>(m, "Gamma")
        .def(py::init<real_t, real_t>(), py::arg("k") = 1.0, py::arg("mu") = 0.1)
        .def_property("k", &mygetter<0>, &mysetter<0>)
        .def_property("mu", &mygetter<1>, &mysetter<1>);

    py::class_<PiecewiseLinear, Proba, shared_ptr<PiecewiseLinear>>(m, "PiecewiseLinear")
        .def(py::init<RealParams const &, real_t>(), py::arg("theta"), py::arg("step") = 1.0);

    py::class_<Cached, Proba, shared_ptr<Cached>>(m, "Cached")
        .def(py::init<std::shared_ptr<Proba> const &, int>(), py::arg("prob"), py::arg("T"))
        .def_readonly("p", &Proba::theta);

    py::class_<Scaled, Proba, shared_ptr<Scaled>>(m, "Scaled")
        .def(py::init<std::shared_ptr<Proba> const &, real_t>(), py::arg("prob"), py::arg("scale") = 1.0);

    py::class_<PDF, Proba, shared_ptr<PDF>>(m, "PDF")
        .def(py::init<std::shared_ptr<Proba> const &>());


    py::class_<Params>(m, "Params")
        .def(py::init<shared_ptr<Proba> const &, shared_ptr<Proba> const &, real_t, real_t, real_t, real_t, real_t, real_t>(),
                "SIB Params class. prob_i and prob_r parameters are defaults.",
                py::arg("prob_i") = *new Uniform(1.0),
                py::arg("prob_r") = *new Exponential(0.1),
                py::arg("pseed") = 0.01,
                py::arg("psus") = 0.5,
                py::arg("fp_rate") = 0.0,
                py::arg("fn_rate") = 0.0,
                py::arg("pautoinf") = 0.0,
                py::arg("learn_rate") = 0.0)

        .def_readwrite("prob_r", &Params::prob_r)
        .def_readwrite("prob_i", &Params::prob_i)
        .def_readwrite("pseed", &Params::pseed)
        .def_readwrite("psus", &Params::psus)
        .def_readwrite("fp_rate", &Params::fp_rate)
        .def_readwrite("fn_rate", &Params::fn_rate)
        .def_readwrite("pautoinf", &Params::pautoinf)
        .def_readwrite("learn_rate", &Params::learn_rate)
        .def("__repr__", &lexical_cast<string, Params>);

    py::class_<FactorGraph>(m, "FactorGraph", "SIB class representing the graphical model of the epidemics")
        .def(py::init<Params const &,
                vector<tuple<int,int,times_t,real_t>>,
                vector<tuple<int,int,times_t>>,
                vector<tuple<int,shared_ptr<Proba>,shared_ptr<Proba>,shared_ptr<Proba>,shared_ptr<Proba>>>
                >(),
                py::arg("params") = Params(shared_ptr<Proba>(new Uniform(1.0)), shared_ptr<Proba>(new Exponential(0.5)), 0.1, 0.45, 0.0, 0.0, 0.0, 0.0),
                py::arg("contacts") = vector<tuple<int,int,times_t,real_t>>(),
                py::arg("observations") = vector<tuple<int,int,times_t>>(),
                py::arg("individuals") = vector<tuple<int,shared_ptr<Proba>,shared_ptr<Proba>,shared_ptr<Proba>,shared_ptr<Proba>>>())
        .def("update", &FactorGraph::iteration,
                py::arg("damping") = 0.0,
                py::arg("learn") = false,
                "perform one iteration")
        .def("loglikelihood", &FactorGraph::loglikelihood, "compute the bethe log-likelihood")
        .def("__repr__", &lexical_cast<string, FactorGraph>)
        .def("append_contact", (void (FactorGraph::*)(int,int,times_t,real_t,real_t)) &FactorGraph::append_contact,
                py::arg("i"),
                py::arg("j"),
                py::arg("t"),
                py::arg("lambdaij"),
                py::arg("lambdaji") = real_t(FactorGraph::DO_NOT_OVERWRITE),
                "appends a new contact from i to j at time t with transmission probabilities lambdaij, lambdaji")
        .def("reset_observations", &FactorGraph::reset_observations,
                py::arg("obs"),
                "resets all observations")
        .def("append_observation", &FactorGraph::append_observation,
                py::arg("i"),
                py::arg("s"),
                py::arg("t"),
                "appends a new observation with state s to node i at time t")
        .def("drop_contacts", &FactorGraph::drop_contacts, "drop contacts at time t (first time)")
        .def("drop_time", &drop_time, "drop time t (first time)")
        .def("drop_sc", &drop_sc,
                py::arg("t"),
                py::arg("maxit_bp") = 1,
                py::arg("tol_bp") = 1e-3,
                py::arg("damping_bp") = 0.0,
                py::arg("maxit_sc") = 20,
                py::arg("tol_sc") = 1e-3,
                py::arg("damping_sc") = 0.1,
                "drop contacts at time t (first time), adjusting fields")

        .def("showmsg", [](FactorGraph & f){f.show_msg(std::cout);}, "show messages for debugging")
        .def_readonly("nodes", &FactorGraph::nodes, "all nodes in this FactorGraph")
        .def_readonly("params", &FactorGraph::params, "parameters");

    py::class_<Node>(m, "Node", "SIB class representing an individual")
        .def("marginal", &get_marginal, "compute marginal probabilities (pS,pI,pR) corresponding to times n.times[1:]")
        .def("marginal_index", &get_marginal_index, "marginal at a given time (excluding time -1)")
        .def_readwrite("ht", &Node::ht, "external prior on ti")
        .def_readwrite("hg", &Node::hg, "external prior on gi")
        .def_readonly("bt", &Node::bt, "belief on ti")
        .def_readonly("bg", &Node::bg, "belief on gi")
        .def_readonly("err", &Node::err_, "error on update")
        .def_readonly("df_i", &Node::df_i, "gradient on prior_i params")
        .def_readonly("df_r", &Node::df_r, "gradient on prior_r params")
        .def_readonly("times", &Node::times, "event times of this node")
        .def_readonly("index", &Node::index, "node index (deprecated, do not use)")
        .def_readonly("prob_i", &Node::prob_i, "probability of infection as function of t-ti")
        .def_readonly("prob_r", &Node::prob_r, "cumulative probability of recovery P(tr>t)")
        .def_readonly("prob_i0", &Node::prob_i0, "probability of infection as function of t-ti for ti=0")
        .def_readonly("prob_r0", &Node::prob_r0, "cumulative probability of recovery P(tr>t) for ti=0");

    m.def("set_num_threads", &omp_set_num_threads, "sets the maximum number of simultaneous cpu threads");
    m.def("version", [](){return VERSION;}, "compiled version of sib");
}
