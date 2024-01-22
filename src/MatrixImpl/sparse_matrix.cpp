#include "sparse_matrix.h"

#include <algorithm>


namespace sparse_matrix {
	Matrix Matrix::from_triplets(index_t cols, index_t rows, std::vector<Triplet> *triplets)
	{
		std::vector<scalar> vals{};
		std::vector<index_t> col_indx{};
		std::vector<index_t> row_ptr{};

		index_t n_nnz = triplets->size();

		std::sort(triplets->begin(), triplets->end(), [=](const Triplet &t1, const Triplet &t2) {
			index_t i1 = t1.col + t1.row * cols;
		index_t i2 = t2.col + t2.row * cols;
		return i1 < i2;
			});


		index_t prev_col = -1;
		index_t prev_row = -1;

		for (index_t i = 0; i < triplets->size(); i++) {
			const auto &t = triplets->at(i);

			if (t.col == prev_col && t.row == prev_row) {
				vals.at(vals.size() - 1) += t.val;
			}
			else if (t.col != prev_col && t.row == prev_row) {
				vals.push_back(t.val);
				col_indx.push_back(t.col);
			}
			else if (t.row != prev_row) {
				vals.push_back(t.val);
				col_indx.push_back(t.col);
				row_ptr.push_back(vals.size() - 1);
			}

			prev_col = t.col;
			prev_row = t.row;
		}

		row_ptr.push_back(vals.size());

		return Matrix{
			.rows = rows,
			.cols = cols,
			.vals = vals,
			.col_indx = col_indx,
			.row_ptr = row_ptr,
		};
	}

	scalar Matrix::get(index_t row, index_t col)
	{
		db_assert(col < cols &&col >= 0);
		db_assert(row < rows &&row >= 0);

		index_t start = row_ptr.at(row);
		index_t end = row_ptr.at(row + 1);

		for (index_t i = start; i < end; i++) {
			if (col_indx.at(i) == col) {
				return vals.at(i);
			}
		}

		return 0;
	}

}

TEST(sparse_matrix_init, {

	std::vector<Triplet> triplets;
	triplets.push_back({ 2, 1, 16.f });
	triplets.push_back({ 3, 2, 11.f });
	triplets.push_back({ 3, 4, 13.f });
	triplets.push_back({ 0, 0, 10.f });
	triplets.push_back({ 0, 3, 12.f });
	triplets.push_back({ 1, 2, 11.f });
	triplets.push_back({ 1, 4, 13.f });

	auto m = sparse_matrix::Matrix::from_triplets(5, 4, &triplets);

	scalar vals[] = {10, 12, 11, 13, 16, 11, 13};
	index_t col_indx[] = {0, 3, 2, 4, 1, 2, 4};
	index_t row_ptr[] = {0, 2, 4, 5, 7};


	for (i32 i = 0; i < m.vals.size(); i++) {
		test_assert(cmp_scalar(m.vals.at(i), vals[i]));
	}

	for (i32 i = 0; i < m.col_indx.size(); i++) {
		test_assert(col_indx[i] == m.col_indx.at(i));
	}


	for (i32 i = 0; i < m.row_ptr.size(); i++) {
		test_assert(row_ptr[i] == m.row_ptr.at(i));
	}

	for (const auto &t : triplets) {
		test_assert(cmp_scalar(m.get(t.row, t.col), t.val));
	}

return {};
	});