#pragma once

#include <set>
#include "geometry.h"
#include "utils.h"



//* functions / types defined in every FEM implementation *//

std::vector<scalar> solve_fem(const Mesh &mesh);
