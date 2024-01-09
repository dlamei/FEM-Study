#pragma once

#include "../utils.h"

namespace linalg {

    typedef u32 index_t;

    struct Matrix3x3 {
        scalar data[9] {0};
    };

    struct Matrix3x2 {
        scalar data[6] {0};
    };

    struct Matrix2x3 {
        scalar data[6] {0};
    };

    struct Vector3 {
        scalar data[3] {0};
    };

    struct IVector3 {
        i32 data[3] {0};
    };


    Matrix3x3 init_mat3x3(
        scalar, scalar, scalar,
        scalar, scalar, scalar,
        scalar, scalar, scalar
    );

    Matrix2x3 init_mat2x3(
        scalar, scalar, scalar,
        scalar, scalar, scalar
    );

    Matrix3x2 init_mat3x2(
        scalar, scalar, scalar,
        scalar, scalar, scalar
    );

    Vector3 init_fvec3(scalar, scalar, scalar);
    IVector3 init_ivec3(i32, i32, i32);


    inline Matrix3x3 clone(const Matrix3x3 &m) { return m; }
    inline Matrix2x3 clone(const Matrix2x3 &m) { return m; }
    inline Matrix3x2 clone(const Matrix3x2 &m) { return m; }
    inline Vector3   clone(const Vector3   &v) { return v; }
    inline IVector3  clone(const IVector3  &v) { return v; }

    // getters
    scalar get(const Matrix3x3 &m, index_t x, index_t y);
    scalar get(const Matrix2x3 &m, index_t x, index_t y);
    scalar get(const Matrix3x2 &m, index_t x, index_t y);
    scalar get(const Vector3   &v, index_t x);
    scalar get(const IVector3  &v, index_t x);

    // setters
    void set(Matrix3x3 *m, index_t x, index_t y, scalar s);
    void set(Matrix2x3 *m, index_t x, index_t y, scalar s);
    void set(Matrix3x2 *m, index_t x, index_t y, scalar s);

    void elem_mul(Matrix3x3 *m, scalar s);
    void elem_mul(Vector3 *v, scalar s);

    Matrix3x3 mat3x2_mul_mat2x3(const Matrix3x2 &m1, const Matrix2x3 &m2);
}
