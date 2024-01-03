#include "FEM.h"


scalar source_fn(const Vector<2> &x) {
    return (8.0 * M_PI * M_PI + 1) 
        * std::cos(2 * M_PI * x(0)) 
        * std::cos(2 * M_PI * x(1));
}


int main() {

    auto mesh = Mesh::load("../meshes/first_mesh_n5.msh");
    auto res = solve_fem(mesh, &source_fn);

    return 0;
}

