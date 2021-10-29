#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ising.hpp"
#include "percolation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(engine, m) {
  py::class_<ising>(m, "ising")
    .def(py::init<uint32_t, uint32_t, double>())
    .def("update", &ising::update)
    .def("config", &ising::config);
  py::class_<percolation>(m, "percolation")
    .def(py::init<uint32_t, uint32_t, double>())
    .def("update", &percolation::update)
    .def("config", &percolation::config);

}
