#include "FEM_naive.h"
#include "benchmark.h"




int main() {

    std::string file_name = "32n_1h";
    std::string folder_name = "1h_exponential_series/";

    // mesh
    auto mesh = Mesh::parse_mesh("../meshes/" + folder_name + file_name + ".msh");

    // Eigen FEM
    auto res = solve_fem(mesh);
    mesh.save_mesh_3D(file_name + "_naive.vtk", res.data, res.height); 


#if PROFILING
    benchmark::global_timer::write_to_file("benchmark.json");
#endif

    return 0;
}

