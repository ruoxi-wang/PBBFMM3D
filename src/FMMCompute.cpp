#include <compute.hpp>
#include "bbfmm.h"
#include "mykernel_python.hpp"
#include <H2_3D_Tree.hpp>
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

using namespace boost::python;


void convert_to_numpy(const vector<double> & input, object obj)
{
    PyObject* pobj = obj.ptr();
    Py_buffer pybuf;
    PyObject_GetBuffer(pobj, &pybuf, PyBUF_SIMPLE);
    void *buf = pybuf.buf;
    double *p = (double*)buf;
    Py_XDECREF(pobj);

    for (int i  = 0; i < input.size(); i++)
    {
        p[i] = input[i];
    }
}

void convert_to_vecOfdouble(const object& obj, vector<double>& output)
{
    PyObject* pobj = obj.ptr();
    Py_buffer pybuf;
    PyObject_GetBuffer(pobj, &pybuf, PyBUF_SIMPLE);
    void *buf = pybuf.buf;
    double *p = (double*)buf;
    Py_XDECREF(pobj);
    output.resize(len(obj));

    for (int i  = 0; i < len(obj); i++)
    {
        output[i] = p[i];
    }
}

void convert_to_vecOfvec3(const object& obj, vector<vector3>& output)
{
    PyObject* pobj = obj.ptr();
    Py_buffer pybuf;
    PyObject_GetBuffer(pobj, &pybuf, PyBUF_SIMPLE);
    void *buf = pybuf.buf;
    double *p = (double*)buf;
    Py_XDECREF(pobj);
    output.resize(len(obj));

    for (int i  = 0; i < len(obj); i++)
    {
        output[i].x = p[i*3+0];
        output[i].y = p[i*3+1];
        output[i].z = p[i*3+2];
    }
}


BOOST_PYTHON_MODULE(FMMCompute)
{
   PyEval_InitThreads();
   class_<std::vector<vector3> >("vecOfvec3")
    .def(vector_indexing_suite<std::vector<vector3> >())
    ;
   class_<std::vector<double> >("vecOfdouble")
    .def(vector_indexing_suite<std::vector<double> >())
    ;
   class_<H2_3D_Compute<myKernel>, boost::noncopyable>("Compute", init<myKernel& , std::vector<vector3>& , std::vector<vector3>&,  std::vector<double>& ,int , std::vector<double>& >())
   .def(init<myKernel& , std::vector<vector3>& , std::vector<vector3>&, std::vector<double>& ,int , std::vector<double>& >())
        ;
   def("convert_to_numpy", convert_to_numpy);
   def("convert_to_vecOfdouble", convert_to_vecOfdouble);
   def("convert_to_vecOfvec3", convert_to_vecOfvec3);
}   

