#pragma once

#include "../utils.h"

/* if this implementation is chosen this will be the matrix type */
typedef struct NaiveMatrix Matrix;

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
struct NaiveMatrix {
	scalar *data{nullptr};
	usize width {0}, height {0};

	// get number of elements
	inline usize count() const;

	// retrieve element
	scalar get(usize indx) const;
	scalar get(usize x, usize y) const;

	// set element
	void set(usize indx, scalar val);
	void set(usize x, usize y, scalar val);

	/* operators */

	// in-place element-wise addition
	static void add_assign(NaiveMatrix *a, const NaiveMatrix &b);
	// add and return result
	static NaiveMatrix add(const NaiveMatrix &a, const NaiveMatrix &b);
	// matrix multiplication
	static NaiveMatrix mul(const NaiveMatrix &a, const NaiveMatrix &b);
	
	static bool eq(const NaiveMatrix &a, const NaiveMatrix &b, float eps = SCALAR_EPS);

	/* constructors / destructors */

	// zero matrix with dim: (w, h)
	static NaiveMatrix zero(usize w, usize h);
	// identity matrix with dim: (w, h)
	static NaiveMatrix ident(usize w, usize h);
	// matrix from a 1D array in column major form
	static NaiveMatrix from_arr(usize w, usize h, const std::initializer_list<scalar> &arr);
	// create matrix from a deep-copy
	static NaiveMatrix clone(const NaiveMatrix &m);

	// deallocate matrix memory
	static void destroy(NaiveMatrix *mat);
};

