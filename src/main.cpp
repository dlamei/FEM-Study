#include "FEM.h"

#include "benchmark.h"
#include "Geometries/geometry_eigen.h"

#define PI 3.14159

scalar source_fn(const Vector<2> &x) {
    return 0;
}


int main() {


    auto mesh = Mesh::load("../meshes/first_mesh_n100.msh");

    std::string file_name = "20n_1h";
    std::string folder_name = "1h_linear_series/";
    auto mesh = Geometry::Mesh::parse_mesh("../meshes/" + folder_name + file_name + ".msh");
    auto res = solve_fem(mesh, &source_fn);


    benchmark::global_timer::write_to_file("benchmark.json");


    mesh.save_mesh_3D(file_name + ".vtk", res);
    return 0;
}

