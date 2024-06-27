#include <pybind11/pybind11.h>
#include <pybind11/eigen/matrix.h>
#include <pybind11/stl.h>
#include <Eigen/Core>
#include <map>
#include <string>

#define STRINGIFY(x) #x

namespace pyeitsolver
{
    bool initialized = false;

    std::map<std::string, int> init(std::string meshfilename, std::string currentfilename)
    {
        initialized = true;
        // Return problem data
        return std::map<std::string, int>{
            {"coeffCount", 1047},
            {"currentsCount", 32},
            {"electrodesCount", 32},
            {"groundNode", 1046}};
    }

    bool is_initialized()
    {
        return initialized;
    }

    std::pair<int, Eigen::VectorXd> solve_forward_problem(const Eigen::VectorXd &conds)
    {
        return std::pair<int, Eigen::VectorXd>(0, conds);
    }

    Eigen::VectorXd solve_full_forward_problem(const Eigen::VectorXd &conds)
    {
        return conds;
    }

}

namespace py = pybind11;

PYBIND11_MODULE(_core, m)
{
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";
    m.def("init", &pyeitsolver::init);
    m.def("is_initialized", &pyeitsolver::is_initialized);
    m.def("solve_forward_problem", &pyeitsolver::solve_forward_problem);
    m.def("solve_full_forward_problem", &pyeitsolver::solve_full_forward_problem);

    #ifdef VERSION_INFO
        m.attr("__version__") = STRINGIFY(VERSION_INFO);
    #else
        m.attr("__version__") = "dev";
    #endif
}