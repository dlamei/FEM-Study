#include "stack_matrix.h"

namespace linalg {

    inline Matrix3x3 init_mat3x3(
            scalar s11, scalar s21, scalar s31,
            scalar s12, scalar s22, scalar s32,
            scalar s13, scalar s23, scalar s33)
    {
        return Matrix3x3 {
            .data = {
                s11, s21, s31,  
                s12, s22, s32, 
                s13, s23, s33 
            },
        };
    }

    inline Matrix2x3 init_mat2x3(
            scalar s11, scalar s21, scalar s31,
            scalar s12, scalar s22, scalar s32) 
    {
        return Matrix2x3 {
            .data = {
                s11, s21, s31,  
                s12, s22, s32, 
            },
        };
    }

    inline Matrix3x2 init_mat3x2(
            scalar s11, scalar s21,
            scalar s12, scalar s22,
            scalar s13, scalar s23)
    {
        return Matrix3x2 {
            .data = {
                s11, s21,
                s12, s22,
                s13, s23
            },
        };

    }

    inline Vector3 init_fvec3(scalar s1, scalar s2, scalar s3) {
        return Vector3 {
            .data = { s1, s2, s3 },
        };
    }

    inline IVector3 init_ivec3(i32 i1, i32 i2, i32 i3) {
        return IVector3 {
            .data = { i1, i2, i3 },
        };
    };

    inline void set(Matrix3x3 *m, index_t x, index_t y, scalar s) {
        assert(x < 3);
        assert(y < 3);

        index_t index = x + y * 3;
        m->data[index] = s;
    };

    inline void set(Matrix2x3 *m, index_t x, index_t y, scalar s) {
        assert(x < 2);
        assert(y < 3);

        index_t index = x + y * 2;
        m->data[index] = s;
    };

    inline void set(Matrix3x2 *m, index_t x, index_t y, scalar s) {
        assert(x < 3);
        assert(y < 2);

        index_t index = x + y * 3;
        m->data[index] = s;
    };

    inline scalar get(const Matrix3x3 &m, index_t x, index_t y) {
        assert(x < 3);
        assert(y < 3);

        index_t index = x + y * 3;
        return m.data[index];
    };

    inline scalar get(const Matrix2x3 &m, index_t x, index_t y) {
        assert(x < 2);
        assert(y < 3);

        index_t index = x + y * 2;
        return m.data[index];
    };

    inline scalar get(const Matrix3x2 &m, index_t x, index_t y) {
        assert(x < 3);
        assert(y < 2);

        index_t index = x + y * 3;
        return m.data[index];
    };

    inline scalar get(const Vector3   &v, index_t x) {
        assert(x < 3);
        return v.data[x];
    };

    inline scalar get(const IVector3  &v, index_t x) {
        assert(x < 3);
        return v.data[x];
    };


    inline void elem_mul(Matrix3x3 *m, scalar s) {
        for (u32 i = 0; i < 9; i++) {
            m->data[i] *= s;
        }
    };

    inline void elem_mul(Vector3 *v, scalar s) {
        v->data[0] *= s;
        v->data[1] *= s;
        v->data[2] *= s;
    };

    inline scalar row_mul_col(const Matrix3x2 &m1, u32 row, const Matrix2x3 &m2, u32 col) {
        return  get(m1, 0, row) * get(m2, col, 0) +
                get(m1, 1, row) * get(m2, col, 1) +
                get(m1, 2, row) * get(m2, col, 2);
    }

    Matrix3x3 mat3x2_mul_mat2x3(const Matrix3x2 &m1, const Matrix2x3 &m2) {
        scalar m11 = row_mul_col(m1, 0, m2, 0);
        scalar m12 = row_mul_col(m1, 1, m2, 0);
        scalar m13 = row_mul_col(m1, 2, m2, 0);

        scalar m21 = row_mul_col(m1, 0, m2, 1);
        scalar m22 = row_mul_col(m1, 1, m2, 1);
        scalar m23 = row_mul_col(m1, 2, m2, 1);

        scalar m31 = row_mul_col(m1, 0, m2, 2);
        scalar m32 = row_mul_col(m1, 1, m2, 2);
        scalar m33 = row_mul_col(m1, 2, m2, 2);

        return init_mat3x3(
                m11, m21, m31,
                m12, m22, m32,
                m13, m23, m33
                );
    };

} // namespace linalg
