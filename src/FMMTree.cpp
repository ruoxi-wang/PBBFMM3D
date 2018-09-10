#include <H2_3D_Tree.hpp>
#include <mykernel_python.hpp>
#include <boost/python.hpp>
using namespace boost::python;


class BaseWrap : public H2_3D_Tree, public wrapper<H2_3D_Tree>
{
    public:
    BaseWrap(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){}

    double EvaluateKernel(vector3& targetpos, vector3& sourcepos)
    {
         return this->get_override("EvaluateKernel")(boost::ref(targetpos), boost::ref(sourcepos));
    }
};

BOOST_PYTHON_MODULE(FMMTree)

{
    PyEval_InitThreads();

    class_<vector3>("vector3", init<double, double, double>())
    .def(init<>())
    .def_readwrite("x", &vector3::x)
    .def_readwrite("y", &vector3::y)
    .def_readwrite("z", &vector3::z);
    ;
    class_<myKernel, boost::noncopyable >("myKernel", init<double, int, int, double, int>())
    .def("buildFMMTree", &H2_3D_Tree::buildFMMTree)
    .def("EvaluateKernel", &H2_3D_Tree::EvaluateKernel)
    ;
   
    class_<BaseWrap, boost::noncopyable>("Atree", init<double, int, int, double, int>())
    .def("EvaluateKernel", pure_virtual(&H2_3D_Tree::EvaluateKernel))
    .def("buildFMMTree", &H2_3D_Tree::buildFMMTree)
    ;
}   
