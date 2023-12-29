#pragma once

#include "utils.h"

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
	scalar *data{NULL};
	usize width {0}, height {0};

	// get number of elements
	inline usize count() const;

	// retrieve element
	scalar get(usize indx) const;
	scalar get(usize x, usize y) const;

	// set element
	void set(scalar val, usize indx);
	void set(scalar val, usize x, usize y);

	/* operators */

	// in-place element-wise addition
	void add_assign(Matrix *a, const Matrix &b);
	// add and return result
	Matrix add(const Matrix &a, const Matrix &b);
	// matrix multiplication
	Matrix mul(const Matrix &a, const Matrix &b);

	/* constructors / destructors */

	// zero matrix with dim: (w, h)
	static Matrix zero(usize w, usize h);
	// identity matrix with dim: (w, h)
	static Matrix ident(usize w, usize h);
	// create matrix from a deep-copy
	static Matrix clone(const Matrix &m);

	// deallocate matrix memory
	static void destroy(Matrix *mat);
};
