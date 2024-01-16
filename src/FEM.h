#pragma  once

#include <Eigen/Eigen>
#include <Eigen/Dense>
#include "geometry.h"
#include "utils.h"


//* functions / types defined in every FEM implementation *//

// store generic Geometry in a implementation specific struct
struct Mesh;
// function pointer for generating load vector
typedef scalar (*source_fn_ptr)(const Eigen::Vector<scalar, 2> &);


Eigen::Vector<scalar, Eigen::Dynamic> solve_fem(const Geometry &mesh, source_fn_ptr source_fn);
