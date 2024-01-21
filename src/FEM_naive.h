#pragma  once

#include <set>
#include "MatrixImpl/naive_matrix.h"
#include "MatrixImpl/stack_matrix.h"
#include "geometry.h"
#include "utils.h"



    //* functions / types defined in every FEM implementation *//

    naive_matrix::Matrix solve_fem(const Mesh &mesh);

