
#pragma once

#include "../utils.h"
#include "../matrix.h"

namespace matrix_impl {

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
    struct Naive {
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
        static void add_assign(Naive *a, const Naive &b);
        // add and return result
        static Naive add(const Naive &a, const Naive &b);
        // matrix multiplication
        static Naive mul(const Naive &a, const Naive &b);

        static bool eq(const Naive &a, const Naive &b, float eps = SCALAR_EPS);

        /* constructors / destructors */

        // zero matrix with dim: (w, h)
        static Naive zero(usize w, usize h);
        // identity matrix with dim: (w, h)
        static Naive ident(usize w, usize h);
        // matrix from a 1D array in column major form
        static Naive from_arr(usize w, usize h, const std::initializer_list<scalar> &arr);
        // create matrix from a deep-copy
        static Naive clone(const Naive &m);

        // deallocate matrix memory
        static void destroy(Naive *mat);
    };
}
