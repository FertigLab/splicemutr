
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
#include <pybind11/stl.h>
#include "splicemutr.hpp"

PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);

namespace py = pybind11;

PYBIND11_MODULE(splicemutrcpp, m){ //for some reason this works. 
    m.doc() = "cpp portion of splicemutr"; // optional module docstring
    m.def("splicemutr", &splicemutr, "running splicemutrcpp");
};