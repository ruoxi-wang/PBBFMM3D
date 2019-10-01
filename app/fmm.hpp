#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include "bbfmm3d.hpp"

class MatrixFMM;
using Eigen::SparseMatrix;
namespace Eigen {
namespace internal {
  // MatrixFMM looks-like a SparseMatrix, so let's inherits its traits:
  template<>
  struct traits<MatrixFMM> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
  {};
}
}
// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixFMM : public Eigen::EigenBase<MatrixFMM> {
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
  Eigen::Product<MatrixFMM,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
    return Eigen::Product<MatrixFMM,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
  }
  // Custom API:
  MatrixFMM() {}
  void initialize(double L_,  int level_, int order_, 
                  double eps_, int use_chebyshev_,
                  std::vector<vector3>& target_, std::vector<vector3>& source_) {
    assert(target_.size() == source_.size());
    this->N = source_.size();
    this->target = target_;
    this->source = source_;
    this->L = L_;
    this->level = level_;
    this->order = order_;
    this->eps = eps_;
    this->use_chebyshev = use_chebyshev_;
  }
public: // for ease of coding
  int N; // problem/matrix size
  std::vector<vector3> target;
  std::vector<vector3> source;
  double L, eps;
  int level, order, use_chebyshev;
};
// Implementation of MatrixFMM * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
namespace internal {
  template<typename Rhs>
  struct generic_product_impl<MatrixFMM, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixFMM,Rhs,generic_product_impl<MatrixFMM,Rhs> >
  {
    typedef typename Product<MatrixFMM,Rhs>::Scalar Scalar;
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const MatrixFMM& lhs, const Rhs& rhs, const Scalar& alpha)
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
      // create tree from scratch
      kernel_Laplacian FMMtree(lhs.L, lhs.level, lhs.order, lhs.eps, lhs.use_chebyshev);
      FMMtree.buildFMMTree();
      H2_3D_Compute<kernel_Laplacian> compute(FMMtree, lhs.target, lhs.source, weight, 1, output_dir);
      for (int i=0; i<N; i++)
          dst(i) = output_dir[i]+1e4*rhs(i); // add diagonal for fast convergence
    }
  };
}
}
