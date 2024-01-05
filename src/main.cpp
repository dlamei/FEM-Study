#include "FEM.h"
#include "Geometries/geometry_eigen.h"

#define PI 3.14159

/*
scalar source_fn(const Vector<2> &x) {
    return (8.0 * PI * PI + 1) 
        * std::cos(2 * PI * x(0)) 
        * std::cos(2 * PI * x(1));
}
*/
scalar source_fn(const Vector<2> &x) {
    return 0;
}

int main() {
    
    std::string file_name = "first_mesh_n5";
    auto mesh = Geometry::Mesh::parse_mesh("../meshes/" + file_name + ".msh");
    auto res = solve_fem(mesh, &source_fn);
    
    mesh.save_mesh_3D(file_name + ".vtk", res);
    return 0;
}

