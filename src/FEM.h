#pragma  once

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include "Geometries/geometry_eigen.h"
#include "utils.h"

namespace Dim {
    const int Dynamic = Eigen::Dynamic;
}

template <const int rows, const int cols>
using Matrix = Eigen::Matrix<scalar, rows, cols>;

template <const int rows>
using Vector = Eigen::Matrix<scalar, rows, 1>;

// QUESTION: When to use "using" and when "typedef"

typedef Matrix<2, 3> Triangle;
typedef Eigen::SparseMatrix<scalar> SparseMatrix;

// QUESTION: WTF is this?
typedef scalar (*source_fn_ptr)(const Vector<2> &);

Vector<Dim::Dynamic> solve_fem(const Geometry::Mesh &mesh, source_fn_ptr source_fn);
