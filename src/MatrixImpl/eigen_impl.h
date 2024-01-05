#pragma once

#include "../utils.h"
#include <Eigen/Eigen>

namespace matrix_impl {

    template <const int rows, const int cols>
    using Matrix = Eigen::Matrix<scalar, rows, cols>;

    template <const int rows>
    using Vector = Eigen::Matrix<scalar, rows, 1>;

    template <const int rows>
    using IVector = Eigen::Matrix<i32, rows, 1>;

    typedef Eigen::SparseMatrix<scalar> SparseMatrix;
}

// map linalg types to eigen

typedef matrix_impl::Matrix<3, 3>   Matrix3x3;
typedef matrix_impl::Matrix<2, 3>   Matrix2x3;
typedef matrix_impl::Matrix<3, 2>   Matrix3x2;
typedef matrix_impl::Vector<3>      Vector3;
typedef matrix_impl::IVector<3>     IVector3;

typedef u32 index_t;
typedef Eigen::Triplet<scalar, index_t> Triplet;

typedef matrix_impl::SparseMatrix SparseMatrix;
typedef matrix_impl::Vector<Eigen::Dynamic> VectorDyn;


// map linalg functions to eigen
namespace linalg {

    inline Matrix3x3 init_mat3x3(
        scalar m11, scalar m21, scalar m31,
        scalar m12, scalar m22, scalar m32,
        scalar m13, scalar m23, scalar m33)
    {
        Matrix3x3 m;
        m << m11, m21, m31,
             m12, m22, m32,
             m13, m23, m33;

        return m;
    }

    inline Matrix2x3 init_mat2x3(
        scalar m11, scalar m21, scalar m31,
        scalar m12, scalar m22, scalar m32)
    {
        Matrix2x3 m;
        m << m11, m21, m31,
             m12, m22, m32;

        return m;
    }

    inline SparseMatrix sparse_from_triplets(u32 w, u32 h, Triplet *begin, usize n) {
        SparseMatrix sparse(w, h);
        sparse.setFromTriplets(begin, begin + n);
        return sparse;
    }

    inline Vector3 init_fvec3(scalar v1, scalar v2, scalar v3) {
        Vector3 v;
        v << v1, v2, v3;
        return v;
    }

    inline IVector3 init_ivec3(i32 v1, i32 v2, i32 v3) {
        IVector3 v;
        v << v1, v2, v3;
        return v;
    }

    inline Matrix3x3 clone(const Matrix3x3 &m) { return m; };
    inline Matrix2x3 clone(const Matrix2x3 &m) { return m; };
    inline Matrix3x2 clone(const Matrix3x2 &m) { return m; };
    inline Vector3   clone(const Vector3   &v) { return v; };
    inline IVector3  clone(const IVector3  &v) { return v; };
    inline VectorDyn clone(const VectorDyn &v) { return v; };

    inline void destroy(SparseMatrix *) {}
    inline void destroy(VectorDyn *) {}

    inline scalar get(const Matrix3x3 &m, index_t x, index_t y) { return m(x, y); };
    inline scalar get(const Matrix2x3 &m, index_t x, index_t y) { return m(x, y); };
    inline scalar get(const Matrix3x2 &m, index_t x, index_t y) { return m(x, y); };

    inline void set(Matrix3x3 *m, index_t x, index_t y, scalar v) { (*m)(x, y) = v; };
    inline void set(Matrix2x3 *m, index_t x, index_t y, scalar v) { (*m)(x, y) = v; };
    inline void set(Matrix3x2 *m, index_t x, index_t y, scalar v) { (*m)(x, y) = v; };

    inline void elem_mul(Matrix3x3 *m, scalar v) { (*m) *= v; };
    inline void elem_mul(Vector3 *vec, scalar v) { (*vec) *= v; };

    inline Matrix3x3 mat3x2_mul_mat2x3(const Matrix3x2 &m1, const Matrix2x3 &m2) {
        return m1 * m2;
    }


    // solve functions
    inline void inverse_inplace(Matrix3x3 *m) {
        auto inv = m->inverse();
        *m = inv;
    }

    inline VectorDyn solve(const SparseMatrix &m, const VectorDyn &v) {
        Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> solver;
        solver.analyzePattern(m);
        solver.factorize(m);
        return solver.solve(v);
    }
}


