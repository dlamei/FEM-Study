#include "stack_matrix.h"

namespace stack_matrix {

    inline Matrix3x3 Matrix3x3::init(
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

    inline scalar Matrix3x3::get(index_t x, index_t y) const {
        db_assert(x < 3);
        db_assert(y < 3);
        return data[x + y * 3];
    }

    inline void Matrix3x3::set(index_t x, index_t y, scalar s) {
        db_assert(x < 3);
        db_assert(y < 3);
        data[x + y * 3] = s;
    }

    inline void Matrix3x3::mul(scalar s) {
        for (u32 i = 0; i < 9; i++) data[i] *= s;
    }


    inline Matrix3x2 Matrix3x2::init(
            scalar s11, scalar s21, scalar s31,
            scalar s12, scalar s22, scalar s32) 
    {
        return Matrix3x2 {
            .data = {
                s11, s21, s31,  
                s12, s22, s32, 
            },
        };
    }

    inline scalar Matrix3x2::get(index_t x, index_t y) const {
        db_assert(x < 3);
        db_assert(y < 2);
        return data[x + y * 3];
    }

    inline void Matrix3x2::set(index_t x, index_t y, scalar s) {
        db_assert(x < 3);
        db_assert(y < 2);
        data[x + y * 3] = s;
    }

    inline Matrix2x3 Matrix2x3::init(
            scalar s11, scalar s21,
            scalar s12, scalar s22,
            scalar s13, scalar s23)
    {
        return Matrix2x3 {
            .data = {
                s11, s21,
                s12, s22,
                s13, s23
            },
        };
    }

    inline scalar Matrix2x3::get(index_t x, index_t y) const {
        db_assert(x < 2);
        db_assert(y < 3);
        return data[x + y * 2];
    }

    inline void Matrix2x3::set(index_t x, index_t y, scalar s) {
        db_assert(x < 2);
        db_assert(y < 3);
        data[x + y * 2] = s;
    }

    inline Matrix3x2 Matrix2x3::transpose() const {
        return Matrix3x2::init(
                get(0, 0), get(0, 1), get(0, 2),
                get(1, 0), get(1, 2), get(1, 2)
                );
    }

    inline scalar row_mul_col(const Matrix3x2 &m1, u32 row, const Matrix2x3 &m2, u32 col) {
        return 
            m1.get(0, row) * m2.get(col, 0) +
            m1.get(0, row) * m2.get(col, 0) +
            m1.get(0, row) * m2.get(col, 0);
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

        return Matrix3x3::init(
                m11, m21, m31,
                m12, m22, m32,
                m13, m23, m33
                );
    };

    inline void Vector3::mul(scalar s) {
        x *= s;
        y *= s;
        z *= s;
    };

} // namespace stack_matrix
