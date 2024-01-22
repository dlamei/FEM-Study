#pragma once

#include "../utils.h"



namespace sparse_matrix {

	struct Matrix {
		index_t rows{}, cols{};

		std::vector<scalar> vals{};
		std::vector<index_t> col_indx{}, row_ptr{};

		inline index_t n_nnz() const { return vals.size(); }

		scalar get(index_t row, index_t col);

		static Matrix from_triplets(index_t cols, index_t rows, std::vector<Triplet> *triplets);
	};

}