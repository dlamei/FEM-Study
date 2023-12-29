#include "linalg.h"

#include <cstring>

/* macros because then the assert then also prints out the function name */

// coordinates to column major index
#define TO_INDEX(indx_var, x, y) \
	assert(x < this->width); \
	assert(y < this->height); \
	usize indx_var = x + y * this->width;

#define CHECK_INDX_BOUND(indx) assert(indx < this->width * this->height);

// DISCLAIMER: only allocates memory for a matrix
Matrix uninit_matrix(usize w, usize h) {
	Matrix mat{};

	usize size = w * h * sizeof(scalar);
	void *data = malloc(size);
	assert(data != NULL);
	mat.data = (scalar *)data;

	return mat;
}

// allocate + initialize to zero
Matrix init_mat(usize w, usize h) {
	Matrix mat{};

	usize size = w * h * sizeof(scalar);
	void *data = malloc(size);
	assert(data != NULL);
	std::memset(data, 0, size);

	mat.data = (scalar *)data;
	return mat;
}

inline usize Matrix::count() const {
	return this->width * this->height; 
}

inline scalar Matrix::get(usize indx) const {
	CHECK_INDX_BOUND(indx);
	return this->data[indx];
}

inline scalar Matrix::get(usize x, usize y) const {
	TO_INDEX(indx, x, y);
	return this->get(indx);
}

inline void Matrix::set(scalar val, usize indx) {
	CHECK_INDX_BOUND(indx);
	this->data[indx] = val;
}

inline void Matrix::set(scalar val, usize x, usize y) {
	TO_INDEX(indx, x, y);
	this->data[indx] = val;
}

void Matrix::add_assign(Matrix *a, const Matrix &b) {
	assert(a->width == b.width);
	assert(a->height == b.height);
	for (usize indx = 0; indx < a->count(); indx++) {
		scalar sum = a->get(indx) + b.get(indx);
		a->set(sum, indx);
	}
}

Matrix Matrix::add(const Matrix &a, const Matrix &b) {
	auto sum = Matrix::clone(a);
	Matrix::add_assign(&sum, b);
	return sum;
}

Matrix Matrix::mul(const Matrix &a, const Matrix &b) {
	assert(a.width == b.height);
	auto res = uninit_matrix(b.width, a.height);

	for (usize i = 0; i < a.height; i++) {
		for (usize j = 0; j < b.width; j++) {
			scalar sum = 0.f;

			for (usize k = 0; k < b.height; k++) {
				sum += a.get(i, k) * b.get(k, j);
			}

			res.set(sum, i, j);
		}
	}
}

inline Matrix Matrix::zero(usize w, usize h) {
	auto mat = init_mat(w, h);
	return mat;
}

inline Matrix Matrix::ident(usize w, usize h) {
	auto mat = Matrix::zero(w, h);

	usize n = std::min(w, h);
	for (usize i = 0; i < n; i++) {
		mat.set(1.f, i);
	}
	return mat;
}

inline Matrix Matrix::clone(const Matrix &source) {
	auto mat = init_mat(source.width, source.height);
	for (usize i = 0; i < mat.count(); i++) {
		mat.set(source.get(i), i);
	}
	return mat;
}

inline void Matrix::destroy(Matrix *m) {
	assert(m->data != NULL);

	if (m-> width > 0 && m->height > 0) {
		delete m->data;
		m->data = NULL;
		m->width = 0;
		m->height = 0;
	}
}


