#include "FEM.h"
#include "benchmark.h"

scalar source_fn(const Eigen::Vector<scalar, 2> &x) {
    return 0;
}


int main() {

    std::string file_name = "128n_1h";
    std::string folder_name = "1h_exponential_series/";

    auto mesh = Mesh::parse_mesh("../meshes/" + folder_name + file_name + ".msh");
    auto res = solve_fem(mesh, &source_fn);

    mesh.save_mesh_3D(file_name + ".vtk", res.data(), res.rows());

#if PROFILING
    benchmark::global_timer::write_to_file("benchmark.json");
#endif

    return 0;
}

