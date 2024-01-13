#include "stack_matrix.h"

#include <iomanip>

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

    inline bool Matrix3x3::eq(const Matrix3x3 &m1, const Matrix3x3 &m2, scalar eps) {
        for (u32 i = 0; i < 9; i++) {
            if (!cmp_scalar(m1.data[i], m2.data[i], eps)) return false;
        }

        return true;
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

    inline scalar det_2x2(scalar a11, scalar a21, scalar a12, scalar a22) {
        return a11 * a22 -  a21 * a12;
    }

    inline scalar Matrix3x3::det() {
        scalar a = data[0], b = data[1], c = data[2],  
               d = data[3], e = data[4], f = data[5], 
               g = data[6], h = data[7], i = data[8];

        return a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
    }

    void Matrix3x3::inverse_inplace() {
        scalar a = data[0], b = data[1], c = data[2],  
               d = data[3], e = data[4], f = data[5], 
               g = data[6], h = data[7], i = data[8];

        scalar det_m = det();

        data[0] =   (e*i - f*h) / det_m;
        data[1] = - (b*i - c*h) / det_m;
        data[2] =   (b*f - c*e) / det_m;
        data[3] = - (d*i - f*g) / det_m;
        data[4] =   (a*i - c*g) / det_m;
        data[5] = - (a*f - c*d) / det_m;
        data[6] =   (d*h - e*g) / det_m;
        data[7] = - (a*h - b*g) / det_m;
        data[8] =   (a*e - b*d) / det_m;
    };

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

    void Matrix3x3::print() {
        for (u32 i = 0; i < 9; i++) {
            std::cout << std::setw(5) << data[i];

            if ((i + 1) % 3) {
                std::cout << ", ";
            } else {
                std::cout << "\n";
            }
        }
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


    TEST(det_matrix_3x3, {
        scalar det;

        det = Matrix3x3::init(
                1, 5, 3,
                2, 4, 7,
                4, 6, 2).det();

        test_assert(cmp_scalar(74.f, det));


        det = Matrix3x3::init(
                -5, -5, -5,
                 3, -1, -2,
                 4,  2,  1).det();

        test_assert(cmp_scalar(-10.f, det));
        return {};
    })

    TEST(inverse_matrix_3x3, {

        Matrix3x3 m;
        Matrix3x3 inv_m;

        m = Matrix3x3::init(
                 1, 2,  1,
                 3, 4,  1,
                -1, 2, -1);

        m.inverse_inplace();


        inv_m = Matrix3x3::init(
                -3.f/4.f,  1.f/2.f, -1.f/4.f,
                 1.f/4.f,      0.f,  1.f/4.f,
                 5.f/4.f, -1.f/2.f, -1.f/4.f
                );

        test_assert(Matrix3x3::eq(m, inv_m));

        m = Matrix3x3::init(
                 0,  1,  4,
                -1, -1,  1,
                 4, -2,  0);

        m.inverse_inplace();


        inv_m = Matrix3x3::init(
                 1.f/14.f, -2.f/7.f, 5.f/28.f,
                  1.f/7.f, -4.f/7.f, -1.f/7.f,
                 3.f/14.f,  1.f/7.f,  1.f/28.f
                );

        test_assert(Matrix3x3::eq(m, inv_m));

        return {};
    })

    } // namespace stack_matrix
