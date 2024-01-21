#pragma once

#include "../utils.h"

namespace stack_matrix {

	typedef u32 index_t;

	struct Matrix3x3 {
		scalar data[9]{ 0 };

		static Matrix3x3 init(
			scalar, scalar, scalar,
			scalar, scalar, scalar,
			scalar, scalar, scalar
		);

		static bool eq(const Matrix3x3 &m1, const Matrix3x3 &m2,
			scalar eps = SCALAR_EPS);

		scalar get(index_t x, index_t y) const;
		void set(index_t x, index_t y, scalar s);
		void mul(scalar s);

		scalar det();
		void inverse_inplace();

		void print();
	};

	struct Matrix3x2 {
		scalar data[6]{ 0 };

		static Matrix3x2 init(
			scalar, scalar, scalar,
			scalar, scalar, scalar
		);

		scalar get(index_t x, index_t y) const;
		void set(index_t x, index_t y, scalar s);

		void print();
	};

	struct Matrix2x3 {
		scalar data[6]{ 0 };

		static Matrix2x3 init(
			scalar, scalar,
			scalar, scalar,
			scalar, scalar
		);

		scalar get(index_t x, index_t y) const;
		void set(index_t x, index_t y, scalar s);

		Matrix3x2 transpose() const;

		void print();
	};

	struct Vector3 {
		scalar x{ 0 }, y{ 0 }, z{ 0 };

		void mul(scalar s);
	};

	struct Vector2 {
		scalar x{ 0 }, y{ 0 };
	};

	struct IVector3 {
		i32 x{ 0 }, y{ 0 }, z{ 0 };
	};

	Matrix3x3 mat3x2_mul_mat2x3(const Matrix3x2 &m1, const Matrix2x3 &m2);
}