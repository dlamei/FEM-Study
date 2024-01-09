
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
        scalar *data{nullptr};
        usize width {0}, height {0};

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
    };
}
