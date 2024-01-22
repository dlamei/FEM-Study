#pragma once

#include <set>
#include "sparse_matrix.h"
#include "../utils.h"

namespace fem_ocl {

	std::vector<scalar> run_conj_grad(sparse_matrix::Matrix &A, std::vector<scalar> &b);

}
