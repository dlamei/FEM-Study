#include "band_matrix.h"


namespace band_matrix {

	BandMatrix BandMatrix::zero(index_t mat_width, index_t band_width)
	{
		assert(band_width != 0, "use vector for matrix with band width 1");
		std::vector<scalar> data(band_width * mat_width);

		return BandMatrix{
			.mat_width = mat_width,
			.band_width = band_width,
			.data = data,
		};
	}

	BandMatrix BandMatrix::from_triplets(index_t width, index_t height, const std::vector<Triplet> &triplts) {

		assert(width == height, "only square matricies can be a band matrix");

		index_t max_diag = 0;

		for (const auto &t : triplts) {
			max_diag = std::max(max_diag, std::abs(t.col - t.row));
		}

		index_t band_width = max_diag * 2 + 1;
		auto bm = BandMatrix::zero(width, band_width);

		for (const auto &t : triplts) {
			// only add the lower half
			if (t.row - t.col >= 0) {
				bm.add(t.col, t.row, t.val);
			}
		}

		return bm;
	}

	inline bool in_band(index_t x, index_t y, index_t band_width) {
		return std::abs(x - y) * 2 + 1 <= band_width;
	}

	inline index_t to_index(index_t x, index_t y, index_t band_width) {
		db_assert(in_band(x, y, band_width));

		if (y - x < 0) {
			// move to lower tri matrix
			std::swap(x, y);
		}

		if (y + 1 >= band_width) {
			x -= y + 2 - band_width;
		}

		index_t data_width = (band_width - 1) / 2 + 1;


		return x + y * band_width;
	}


	void BandMatrix::set(index_t x, index_t y, scalar v)
	{
		index_t indx = to_index(x, y, band_width);
		data.at(indx) = v;
	}

	scalar BandMatrix::get(index_t x, index_t y)
	{
		if (in_band(x, y, band_width)) {
			return data.at(to_index(x, y, band_width));
		}
		else {
			return 0;
		}
	}

	void BandMatrix::add(index_t x, index_t y, scalar v)
	{
		db_assert(in_band(x, y, band_width));
		data.at(to_index(x, y, band_width)) += v;
	}
}

void test_band() {
	std::vector<Triplet> triplets;

	for (i32 i = 0; i < 10; i++) {
		triplets.push_back({ i, i, (float)i / 10 });
	}

	triplets.push_back({ 5, 6, 3.2 });
	triplets.push_back({ 6, 5, 3.2 });

	auto bm = band_matrix::BandMatrix::from_triplets(10, 10, triplets);


}


TEST(band_matrix, {

	std::vector<Triplet> triplets;

	for (i32 i = 0; i < 10; i++) {
		triplets.push_back({ i, i, (float)i / 10 });
	}

	triplets.push_back({ 5, 6, 3.2 });
	triplets.push_back({ 6, 5, 3.2 });

	auto bm = band_matrix::BandMatrix::from_triplets(10, 10, triplets);


	for (const auto &t : triplets) {
		test_assert(cmp_scalar(bm.get(t.col, t.row), t.val));
	}


	return {};
	});
