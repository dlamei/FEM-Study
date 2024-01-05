#pragma once


// choose a backend 

#define BACKEND_HEADER "MatrixImpl/eigen_impl.h"
#define BACKEND_SRC    "MatrixImpl/eigen_impl.cpp"



#include BACKEND_HEADER


/* TYPES DEFINED BY BACKEND */


// fixed size

static_assert(is_defined<Matrix3x3>);
static_assert(is_defined<Matrix2x3>);
static_assert(is_defined<Matrix3x2>);

static_assert(is_defined<Vector3>);
static_assert(is_defined<IVector3>);

static_assert(is_defined<Triplet>);
static_assert(is_defined<index_t>);


// dynamic size

static_assert(is_defined<SparseMatrix>);
static_assert(is_defined<VectorDyn>);


/* FUNCTIONS DEFINED BY BACKEND */

namespace linalg {


    // constructors
    Matrix3x3 init_mat3x3(
        scalar, scalar, scalar,
        scalar, scalar, scalar,
        scalar, scalar, scalar
    );

    Matrix2x3 init_mat2x3(
        scalar, scalar, scalar,
        scalar, scalar, scalar
    );

    // create sparse matrix from n triplets
    SparseMatrix sparse_from_triplets(u32 w, u32 h, Triplet *, usize n);

    Vector3 init_fvec3(scalar, scalar, scalar);
    IVector3 init_ivec3(i32, i32, i32);

    Matrix3x3 clone(const Matrix3x3 &);
    Matrix2x3 clone(const Matrix2x3 &);
    Matrix3x2 clone(const Matrix3x2 &);
    Vector3   clone(const Vector3   &);
    IVector3  clone(const IVector3  &);
    VectorDyn clone(const VectorDyn &);

    //destructors
    // only dynamically sized elements need to be deleted
    void destroy(SparseMatrix *);
    void destroy(VectorDyn *);


    // getters
    scalar get(const Matrix3x3 &, index_t, index_t);
    scalar get(const Matrix2x3 &, index_t, index_t);
    scalar get(const Matrix3x2 &, index_t, index_t);
    scalar get(const Vector3   &, index_t);
    scalar get(const IVector3  &, index_t);
    scalar get(const VectorDyn &, index_t);

    // setters
    void set(Matrix3x3 *, index_t, index_t, scalar);
    void set(Matrix2x3 *, index_t, index_t, scalar);
    void set(Matrix3x2 *, index_t, index_t, scalar);
    void get(Vector3   *, index_t, scalar);
    void get(IVector3  *, index_t, scalar);
    void get(VectorDyn *, index_t, scalar);


    // multiplications
    void elem_mul(Matrix3x3 *, scalar);
    void elem_mul(Vector3 *, scalar);

    Matrix3x3 mat3x2_mul_mat2x3(const Matrix3x2 &, const Matrix2x3 &);


    // solve functions
    void inverse_inplace(Matrix3x3 *);

    VectorDyn solve(const SparseMatrix &, const VectorDyn &);
}

