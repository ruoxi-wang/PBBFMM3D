#include <H2_3D_Tree.hpp>
class myKernel: public H2_3D_Tree {
public:
    myKernel(double L, int level, int n, double epsilon, int use_chebyshev):H2_3D_Tree(L,level,n, epsilon, use_chebyshev){};
    void setKernelProperty() {
        this->homogen = 0;
        this->symmetry = 1;
        this->kernelType = "myKernel";
    } 
    double EvaluateKernel(vector3& fieldpos, vector3& sourcepos) {
        vector3 diff;        
        // Compute exp(-r)
        diff.x = sourcepos.x - fieldpos.x;
        diff.y = sourcepos.y - fieldpos.y;
        diff.z = sourcepos.z - fieldpos.z;
        double r = sqrt(diff.x*diff.x/100. + diff.y*diff.y/100. + diff.z*diff.z/2.25);
        return 0.75*exp(-r);
    };
};
