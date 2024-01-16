#include "naive_matrix.h"

#include <cstring>
#include <iomanip>

#include <algorithm>

using namespace naive_matrix;


/* macros because then the assert also prints out the function name */

// coordinates to column major index
#define TO_INDEX(indx_var, x, y) \
    db_assert(x < this->width);  \
    db_assert(y < this->height); \
    u64 indx_var = x + y * this->width;

#define CHECK_INDX_BOUND(indx) \
    db_assert(indx < this->width * this->height);

// DISCLAIMER: only allocates memory for a matrix
Matrix uninit_matrix(u64 w, u64 h) {
    Matrix mat{};

    u64 size = w * h * sizeof(scalar);
    void *data = malloc(size);
    db_assert(data != NULL);
    mat.width = w;
    mat.height = h;
    mat.data = (scalar *)data;

    return mat;
}

// allocate + initialize to zero
Matrix init_mat(u64 w, u64 h) {
    Matrix mat = uninit_matrix(w, h);
    std::memset(mat.data, 0, w * h * sizeof(scalar));
    return mat;
}


inline u64 Matrix::count() const {
    return this->width * this->height; 
}

inline scalar Matrix::get(u64 indx) const {
    CHECK_INDX_BOUND(indx);
    return this->data[indx];
}

inline scalar Matrix::get(u64 x, u64 y) const {
    TO_INDEX(indx, x, y);
    return this->get(indx);
}

inline void Matrix::set(u64 indx, scalar val) {
    CHECK_INDX_BOUND(indx);
    this->data[indx] = val;
}

inline void Matrix::set(u64 x, u64 y, scalar val) {
    TO_INDEX(indx, x, y);
    this->data[indx] = val;
}

void Matrix::add_assign(Matrix *a, const Matrix &b) {
    db_assert(a->width == b.width);
    db_assert(a->height == b.height);
    for (u64 indx = 0; indx < a->count(); indx++) {
        scalar sum = a->get(indx) + b.get(indx);
        a->set(indx, sum);
    }
}

Matrix Matrix::add(const Matrix &a, const Matrix &b) {
    auto sum = Matrix::clone(a);
    Matrix::add_assign(&sum, b);
    return sum;
}

Matrix Matrix::mul(const Matrix &a, const Matrix &b) {
    db_assert(a.width == b.height);
    auto res = uninit_matrix(b.width, a.height);

    for (u64 i = 0; i < b.width; i++) {
        for (u64 j = 0; j < a.height; j++) {
            scalar sum = 0;
            for (u64 k = 0; k < a.width; k++) {
                sum += a.get(k, j) * b.get(i, k);
            }
            res.set(i, j, sum);
        }
    }

    return res;
}


void Matrix::solve(Matrix *a, Matrix *b) {
 // TODO: implement
}


bool Matrix::eq(const Matrix &a, const Matrix &b, float eps) {
    if (a.width != b.width || a.height != b.height) return false;

    for (u64 i = 0; i < a.count(); i++) {
        if (!cmp_scalar(a.get(i), b.get(i), eps)) return false;
    }

    return true;
}

inline Matrix Matrix::zero(u64 w, u64 h) {
    auto mat = init_mat(w, h);
    return mat;
}

inline Matrix Matrix::ident(u64 w, u64 h) {
    auto mat = Matrix::zero(w, h);

    u64 n = std::min(w, h);
    for (u64 i = 0; i < n; i++) {
        mat.set(i, i, 1.f);
    }
    return mat;
}

inline Matrix Matrix::from_arr(u64 w, u64 h, const std::initializer_list<scalar> &arr) {
    db_assert(w * h == arr.size());
    auto mat = uninit_matrix(w, h);

    u64 indx = 0;
    for (scalar val : arr) {
        mat.set(indx, val);
        indx++;
    }

    return mat;
}

inline Matrix Matrix::clone(const Matrix &source) {
    auto mat = init_mat(source.width, source.height);
    for (u64 i = 0; i < mat.count(); i++) {
        mat.set(i, source.get(i));
    }
    return mat;
}

inline void Matrix::destroy(Matrix *m) {
    db_assert(m->data != nullptr);

    if (m-> width > 0 && m->height > 0) {
        delete m->data;
        m->data = nullptr;
        m->width = 0;
        m->height = 0;
    }
}

void Matrix::print() {
    for (u64 i = 0; i < width; i++) {
        for (u64 j = 0; j < width; j++) {
            std::cout << std::setw(5) << get(i, j) << ", ";
        }
        std::cout << "\n";
    }
}


///////////////// TESTS /////////////////

TEST(matrix_eq, {
    auto m1 = Matrix::zero(120, 312);

    for (u64 indx = 0; indx < m1.count(); indx++) {
        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        m1.set(indx, r);
    }

    auto m2 = Matrix::clone(m1);
    test_assert(Matrix::eq(m1, m2));

    return {};
})


TEST(matrix_mul, {
        {
        auto m1 = Matrix::from_arr(3, 3, {6.f, 2.f, 4.f, -1.f, 4.f, 3.f, -2.f, 9.f, 3.f,});
        auto m2 = Matrix::from_arr(1, 3, {4.f, -2.f, -1.f});
        auto m3 = Matrix::from_arr(1, 3, {16, -15, -29});
        defer(Matrix::destroy(&m1));
        defer(Matrix::destroy(&m2));
        defer(Matrix::destroy(&m3));

        test_assert(Matrix::eq(m3, Matrix::mul(m1, m2)), "mul 1");
        }

        {
        auto m1 = Matrix::from_arr(3, 3, {6.f, 2.f, 4.f, -1.f, 4.f, 3.f, -2.f, 9.f, 3.f,});
        auto m2 = Matrix::from_arr(1, 3, {4.f, -2.f, -1.f});
        auto m3 = Matrix::from_arr(1, 3, {17, -15, -29});
        defer(Matrix::destroy(&m1));
        defer(Matrix::destroy(&m2));
        defer(Matrix::destroy(&m3));

        test_assert(!Matrix::eq(m3, Matrix::mul(m1, m2)), "mul 2");
        }

        return {};
})

TEST(matrix_ident, {
    auto ident = Matrix::ident(5, 5);

    for (u64 x = 0; x < ident.width; x++) {
        for (u64 y = 0; y < ident.width; y++) {

            if (x == y) test_assert(cmp_scalar(1.f, ident.get(x, y)));
            else        test_assert(cmp_scalar(0.f, ident.get(x, y)));
        }
    }

    return {};
})

TEST(matrix_solve, {


    auto a = Matrix::zero(3, 3);
    auto b = Matrix::zero(3, 1);

    a.set(0, 0,  4);
    a.set(1, 0, -4);
    a.set(2, 0,  8);

    a.set(0, 1,  8);
    a.set(1, 1,  4);
    a.set(2, 1, -4);

    a.set(0, 2,  12);
    a.set(1, 2,  -8);
    a.set(2, 2, -12);

    b.set(0, 20);
    b.set(1, 4);
    b.set(2, -40);


    Matrix::solve(&a, &b);

    b.print();


    return {};
})
