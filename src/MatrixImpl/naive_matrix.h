
#pragma once

#include "../utils.h"

namespace naive_matrix {

	/*
	 * f32 column major matrix
	 *
	 * types of indexing:
	 *  - index: n-th element in the column major matrix
	 *  - (x, y): coordinates in the matrix
	 *  	x -> horizontal, y -> vertical
	 *  	(0, 0): top left, (width - 1, height - 1): bottom right
	 *
	 */
	struct Matrix {
		scalar *data{ nullptr };
		u64 width{ 0 }, height{ 0 };

		// get number of elements
		inline usize count() const;

		// retrieve element
		scalar get(usize x, usize y) const;
		scalar get(usize indx) const;

		// set element
		void set(usize indx, scalar val);
		void set(usize x, usize y, scalar val);

		/* operators */

		// in-place element-wise addition
		static void add_assign(Matrix *a, const Matrix &b);
		// add and return result
		static Matrix add(const Matrix &a, const Matrix &b);
		// matrix multiplication
		static Matrix mul(const Matrix &a, const Matrix &b);

		// Precondition: L is lower triangular n*n and b is a n*1 matrix;
		static Matrix forward_substitution(const Matrix &L, const Matrix &b);

		// Precondition: U is upper triangular n*n and b is a n*1 matrix;
		static Matrix backward_substitution(const Matrix &U, const Matrix &b);

		// Precondition: A is a LDL decompose Matrix n*n and b is 1*n
		static void inplace_LDL_decomposition(Matrix *a);

		// Precondition: A is a LDL decompose Matrix n*n and b is 1*n
		static void LDLT_solve(Matrix *a, Matrix *b);

		static void solve(Matrix *a, Matrix *b);

		static bool eq(const Matrix &a, const Matrix &b, float eps = SCALAR_EPS);

		/* constructors / destructors */

		// zero matrix with dim: (w, h)
		static Matrix zero(usize w, usize h);
		// identity matrix with dim: (w, h)
		static Matrix ident(usize w, usize h);
		// matrix from a 1D array in column major form
		static Matrix from_arr(usize w, usize h, const std::initializer_list<scalar> &arr);
		// create matrix from a deep-copy
		static Matrix clone(const Matrix &m);

		// deallocate matrix memory
		static void destroy(Matrix *mat);

		void print();
	};
}