#define FEM_CONJ_GRAD_IMPL 1
#define FEM_EIGEN_IMPL 0
#define FEM_NAIVE_IMPL 0

#if FEM_NAIVE_IMPL
#include "FEM_naive.h"
#elif FEM_EIGEN_IMPL
#include "FEM.h"
#elif FEM_CONJ_GRAD_IMPL
#include "FEM_conj_grad.h"
#endif

#include "benchmark.h"



int main() {

	std::string file_name = "100n_1h";
	std::string folder_name = "1h_linear_series/";

	// mesh
	auto mesh = Mesh::parse_mesh("../meshes/" + folder_name + file_name + ".msh");

	auto start_time = std::chrono::high_resolution_clock::now();

#if FEM_NAIVE_IMPL


#elif FEM_EIGEN_IMPL

	// Eigen FEM
	auto res = solve_fem(mesh);
	mesh.save_mesh_3D(file_name + ".vtk", res.data(), res.rows());


#elif FEM_CONJ_GRAD_IMPL

	auto res = solve_fem(mesh);
	mesh.save_mesh_3D(file_name + ".vtk", res.data(), res.size());

#endif

	auto end_time = std::chrono::high_resolution_clock::now();
	std::cout << end_time - start_time << "\n";

#if PROFILING
	benchmark::global_timer::write_to_file("benchmark.json");
#endif

	return 0;
}

