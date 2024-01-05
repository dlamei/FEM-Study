#include "naive_matrix.h"
#include "../utils.h"

#include <cstring>

/* macros because then the assert also prints out the function name */

// coordinates to column major index
#define TO_INDEX(indx_var, x, y) \
	db_assert(x < this->width); \
	db_assert(y < this->height); \
	usize indx_var = x + y * this->width;

#define CHECK_INDX_BOUND(indx) db_assert(indx < this->width * this->height);

// DISCLAIMER: only allocates memory for a matrix
NaiveMatrix uninit_matrix(usize w, usize h) {
	NaiveMatrix mat{};

	usize size = w * h * sizeof(scalar);
	void *data = malloc(size);
	db_assert(data != NULL);
	mat.width = w;
	mat.height = h;
	mat.data = (scalar *)data;

	return mat;
}

// allocate + initialize to zero
NaiveMatrix init_mat(usize w, usize h) {
	NaiveMatrix mat = uninit_matrix(w, h);
	std::memset(mat.data, 0, w * h * sizeof(scalar));
	return mat;
}

inline usize NaiveMatrix::count() const {
	return this->width * this->height; 
}

inline scalar NaiveMatrix::get(usize indx) const {
	CHECK_INDX_BOUND(indx);
	return this->data[indx];
}

inline scalar NaiveMatrix::get(usize x, usize y) const {
	TO_INDEX(indx, x, y);
	return this->get(indx);
}

inline void NaiveMatrix::set(usize indx, scalar val) {
	CHECK_INDX_BOUND(indx);
	this->data[indx] = val;
}

inline void NaiveMatrix::set(usize x, usize y, scalar val) {
	TO_INDEX(indx, x, y);
	this->data[indx] = val;
}

void NaiveMatrix::add_assign(NaiveMatrix *a, const NaiveMatrix &b) {
	db_assert(a->width == b.width);
	db_assert(a->height == b.height);
	for (usize indx = 0; indx < a->count(); indx++) {
		scalar sum = a->get(indx) + b.get(indx);
		a->set(indx, sum);
	}
}

NaiveMatrix NaiveMatrix::add(const NaiveMatrix &a, const NaiveMatrix &b) {
	auto sum = NaiveMatrix::clone(a);
	NaiveMatrix::add_assign(&sum, b);
	return sum;
}

NaiveMatrix NaiveMatrix::mul(const NaiveMatrix &a, const NaiveMatrix &b) {
	db_assert(a.width == b.height);
	auto res = uninit_matrix(b.width, a.height);

    for (usize i = 0; i < b.width; i++) {
        for (usize j = 0; j < a.height; j++) {
            scalar sum = 0;
            for (usize k = 0; k < a.width; k++) {
                sum += a.get(k, j) * b.get(i, k);
            }
            res.set(i, j, sum);
        }
    }

	return res;
}

bool NaiveMatrix::eq(const NaiveMatrix &a, const NaiveMatrix &b, float eps) {
	if (a.width != b.width || a.height != b.height) return false;

	for (usize i = 0; i < a.count(); i++) {
		if (!cmp_scalar(a.get(i), b.get(i), eps)) return false;
	}

	return true;
}

inline NaiveMatrix NaiveMatrix::zero(usize w, usize h) {
	auto mat = init_mat(w, h);
	return mat;
}

inline NaiveMatrix NaiveMatrix::ident(usize w, usize h) {
	auto mat = NaiveMatrix::zero(w, h);

	usize n = std::min(w, h);
	for (usize i = 0; i < n; i++) {
		mat.set(i, i, 1.f);
	}
	return mat;
}

inline NaiveMatrix NaiveMatrix::from_arr(usize w, usize h, const std::initializer_list<scalar> &arr) {
	db_assert(w * h == arr.size());
	auto mat = uninit_matrix(w, h);

	usize indx = 0;
	for (scalar val : arr) {
		mat.set(indx, val);
		indx++;
	}

	return mat;
}

inline NaiveMatrix NaiveMatrix::clone(const NaiveMatrix &source) {
	auto mat = init_mat(source.width, source.height);
	for (usize i = 0; i < mat.count(); i++) {
		mat.set(i, source.get(i));
	}
	return mat;
}

inline void NaiveMatrix::destroy(NaiveMatrix *m) {
	db_assert(m->data != nullptr);

	if (m-> width > 0 && m->height > 0) {
		delete m->data;
		m->data = nullptr;
		m->width = 0;
		m->height = 0;
	}
}

TEST(matrix_eq, {
    auto m1 = NaiveMatrix::zero(120, 312);
    for (usize indx = 0; indx < m1.count(); indx++) {
        float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        m1.set(indx, r);
    }
    auto m2 = NaiveMatrix::clone(m1);

    test_assert(NaiveMatrix::eq(m1, m2));

    return {};
})

TEST(matrix_mul, {
    {
        auto m1 = NaiveMatrix::from_arr(3, 3, {6.f, 2.f, 4.f, -1.f, 4.f, 3.f, -2.f, 9.f, 3.f,});
        auto m2 = NaiveMatrix::from_arr(1, 3, {4.f, -2.f, -1.f});
        auto m3 = NaiveMatrix::from_arr(1, 3, {16, -15, -29});
        defer(NaiveMatrix::destroy(&m1));
        defer(NaiveMatrix::destroy(&m2));
        defer(NaiveMatrix::destroy(&m3));

        test_assert(NaiveMatrix::eq(m3, NaiveMatrix::mul(m1, m2)), "mul 1");
    }

    {
        auto m1 = NaiveMatrix::from_arr(3, 3, {6.f, 2.f, 4.f, -1.f, 4.f, 3.f, -2.f, 9.f, 3.f,});
        auto m2 = NaiveMatrix::from_arr(1, 3, {4.f, -2.f, -1.f});
        auto m3 = NaiveMatrix::from_arr(1, 3, {17, -15, -29});
        defer(NaiveMatrix::destroy(&m1));
        defer(NaiveMatrix::destroy(&m2));
        defer(NaiveMatrix::destroy(&m3));

        test_assert(!NaiveMatrix::eq(m3, NaiveMatrix::mul(m1, m2)), "mul 2");
    }

    return {};
})

TEST(matrix_ident, {
    auto ident = NaiveMatrix::ident(5, 5);
    for (usize x = 0; x < ident.width; x++) {
        for (usize y = 0; y < ident.width; y++) {
            if (x == y) {
                test_assert(cmp_scalar(1.f, ident.get(x, y)));
            } else {
                test_assert(cmp_scalar(0.f, ident.get(x, y)));
            }
        }
    }
    return {};
})


//static __register_test matrix_test([]() { test_matrix_mul(); });
