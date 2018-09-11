#include <H2_3D_Tree.hpp>
class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int tree_level, int interpolation_order, double epsilon, int use_chebyshev):H2_3D_Tree(L,tree_level,interpolation_order, epsilon, use_chebyshev){};
    void SetKernelProperty() {
        this->homogen = 0;
        this->symmetry = 1;
        this->kernelType = "myKernel";
    } 
    double EvaluateKernel(vector3& targetpos, vector3& sourcepos) {
        vector3 diff;        
        // Compute exp(-r)
        diff.x = sourcepos.x - targetpos.x;
        diff.y = sourcepos.y - targetpos.y;
        diff.z = sourcepos.z - targetpos.z;
        double r = sqrt(diff.x*diff.x/100. + diff.y*diff.y/100. + diff.z*diff.z/2.25);
        return 0.75*exp(-r);
    };
};
