#include "utils.h"

/* generic matrix functions */

// choose matrix implementation
// should be equal to implementation chosen in matrix.cpp
#include "MatrixImpl/naive_matrix.h"

namespace matrix {

	Matrix zero(usize w, usize h);
	Matrix ident(usize w, usize h);

	void set(Matrix *m, usize x, usize y, scalar val);
	Matrix add(const Matrix &m1, const Matrix &m2);
	Matrix mul(const Matrix &m1, const Matrix &m2);
	bool eq(const Matrix &m1, const Matrix &m2);

	void destroy(Matrix *m);

	void print(const Matrix &m1);
}
