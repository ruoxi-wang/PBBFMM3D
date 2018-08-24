#include <compute.hpp>
#include "bbfmm.h"
#include "mykernel_python.hpp"
#include <H2_3D_Tree.hpp>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;


BOOST_PYTHON_MODULE(FMMCompute)
{
   PyEval_InitThreads();
   class_<std::vector<vector3> >("vector_vector3")
    .def(vector_indexing_suite<std::vector<vector3> >())
    ;
   class_<std::vector<double> >("vector_double")
    .def(vector_indexing_suite<std::vector<double> >())
    ;
   class_<H2_3D_Compute<myKernel>, boost::noncopyable>("Compute", init<myKernel& , std::vector<vector3>& , std::vector<vector3>&, int , int , std::vector<double>& ,int , std::vector<double>& >())
   .def(init<myKernel& , std::vector<vector3>& , std::vector<vector3>&, int , int , std::vector<double>& ,int , std::vector<double>& >())
        ;
}   