#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "bbfmm3d.hpp"

class MatrixDirect;
using Eigen::SparseMatrix;
namespace Eigen {
namespace internal {
  // MatrixDirect looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<MatrixDirect> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};
}
}
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixDirect : public Eigen::EigenBase<MatrixDirect> {
public:
  // Required typedefs, constants, and method:
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };
  Index rows() const { return N; }
  Index cols() const { return N; }
  template<typename Rhs>
  Eigen::Product<MatrixDirect,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixDirect,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
  // Custom API:
  MatrixDirect() {}
  void initialize(kernel_Laplacian* FMMtree_, 
                  std::vector<vector3>& target_, std::vector<vector3>& source_) {
    assert(target_.size() == source_.size());
    this->N = source_.size();
    this->target = target_;
    this->source = source_;
    this->FMMtree = FMMtree_;
  }
public: // for ease of coding
  int N; // problem/matrix size
  std::vector<vector3> target;
  std::vector<vector3> source;
  kernel_Laplacian* FMMtree;
};
// Implementation of MatrixDirect * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
  template<typename Rhs>
  struct generic_product_impl<MatrixDirect, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixDirect,Rhs,generic_product_impl<MatrixDirect,Rhs> >
  {
    typedef typename Product<MatrixDirect,Rhs>::Scalar Scalar;
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const MatrixDirect& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      assert(alpha==Scalar(1) && "scaling is not implemented");
      EIGEN_ONLY_USED_FOR_DEBUG(alpha);
      // initialize weights
      int N = lhs.N;
      std::vector<double> weight(N);
      for (int i=0; i<N; i++)
          weight[i] = rhs(i);
      std::vector<double> output_dir(N);
      DirectCalc3D(lhs.FMMtree, lhs.target, lhs.source, weight, 1, output_dir, N);
      for (int i=0; i<N; i++)
          dst(i) = output_dir[i]+1e4*rhs(i); // add diagonal for fast convergence
    }
  };
}
}
