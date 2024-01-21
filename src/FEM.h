#pragma  once

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <set>
#include "geometry.h"
#include "utils.h"


//* functions / types defined in every FEM implementation *//
// function pointer for generating load vector
typedef scalar (*source_fn_ptr)(const Eigen::Vector<scalar, 2> &);

Eigen::Vector<scalar, Eigen::Dynamic> solve_fem(const Mesh &mesh, source_fn_ptr source_fn);

