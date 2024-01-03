#pragma  once

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include "utils.h"

namespace Dim {
    const int Dynamic = Eigen::Dynamic;
}

template <const int rows, const int cols>
using Matrix = Eigen::Matrix<scalar, rows, cols>;

template <const int rows>
using Vector = Eigen::Matrix<scalar, rows, 1>;

typedef Matrix<2, 3> Triangle;
typedef Eigen::SparseMatrix<scalar> SparseMatrix;

typedef scalar (*source_fn_ptr)(const Vector<2> &);


struct Mesh {

    usize n_nodes { 0 };
    usize n_triangles { 0 };

    // 2 scalar per node
    Matrix<Dim::Dynamic, 2> nodes{};

    // 3 vert_indx's per triangle
    Eigen::Matrix<int, Eigen::Dynamic, 3> triangles{};

    static Mesh load(const std::string &file_name);
};


Vector<Dim::Dynamic> solve_fem(const Mesh &mesh, source_fn_ptr source_fn);
