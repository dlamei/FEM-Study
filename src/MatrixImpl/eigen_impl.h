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
        scalar s11, scalar s21, scalar s31,
        scalar s12, scalar s22, scalar s32,
        scalar s13, scalar s23, scalar s33)
    {
        Matrix3x3 m;
        m << s11, s21, s31,
             s12, s22, s32,
             s13, s23, s33;

        return m;
    }

    inline Matrix2x3 init_mat2x3(
        scalar s11, scalar s21, scalar s31,
        scalar s12, scalar s22, scalar s32)
    {
        Matrix2x3 m;
        m << s11, s21, s31,
             s12, s22, s32;

        return m;
    }

    inline SparseMatrix sparse_from_triplets(u32 w, u32 h, Triplet *begin, usize n) {
        SparseMatrix sparse(w, h);
        sparse.setFromTriplets(begin, begin + n);
        return sparse;
    }

    inline Vector3 init_fvec3(scalar s1, scalar s2, scalar s3) {
        Vector3 v;
        v << s1, s2, s3;
        return v;
    }

    inline IVector3 init_ivec3(i32 s1, i32 s2, i32 s3) {
        IVector3 v;
        v << s1, s2, s3;
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
    inline scalar get(const Vector3   &v, index_t x) { return v(x); };
    inline scalar get(const IVector3  &v, index_t x) { return v(x); };
    inline scalar get(const VectorDyn &v, index_t x) { return v(x); };

    inline void set(Matrix3x3 *m, index_t x, index_t y, scalar s) { (*m)(x, y) = s; };
    inline void set(Matrix2x3 *m, index_t x, index_t y, scalar s) { (*m)(x, y) = s; };
    inline void set(Matrix3x2 *m, index_t x, index_t y, scalar s) { (*m)(x, y) = s; };
    inline void set(Vector3   *v, index_t x, scalar s) { (*v)(x) = s; };
    inline void set(IVector3  *v, index_t x, scalar s) { (*v)(x) = s; };
    inline void set(VectorDyn *v, index_t x, scalar s) { (*v)(x) = s; };

    inline void elem_mul(Matrix3x3 *m, scalar s) { (*m) *= s; };
    inline void elem_mul(Vector3   *v, scalar s) { (*v) *= s; };

    inline Matrix3x3 mat3x2_mul_mat2x3(const Matrix3x2 &m1, const Matrix2x3 &m2) { return m1 * m2; }


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


