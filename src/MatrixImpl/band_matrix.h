#pragma once

#include "../utils.h"

namespace band_matrix {

	//struct Triplet {
	//	index_t x{}, y{};
	//	scalar v{};
	//};

	// requires the matrix to be symmetric
	struct BandMatrix {

		index_t mat_width{};
		index_t band_width{};
		std::vector<scalar> data{};


		static BandMatrix zero(index_t mat_width, index_t bandwidth);
		static BandMatrix from_triplets(index_t width, index_t height, const std::vector<Triplet> &triplts);

		void set(index_t x, index_t y, scalar v);
		scalar get(index_t x, index_t y);
		void add(index_t x, index_t y, scalar v);

		inline scalar *storage_ptr() { return data.data(); }
	};


}