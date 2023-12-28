#include <iostream>
#include "utils.h"
#include "./Geometries/geometry.h"

void print(int i) {
    std::cout << i << std::endl;
}
template<typename T>
void print_mat(std::vector<std::vector<T>> mat) {
    for(int i = 0; i < mat.size(); ++i) {
        for(int j = 0; j < mat.at(0).size(); ++j) {
            std::cout << mat.at(i).at(j) << " ";
        }
        std::cout << "\n";
    }
}

int main() {

    //std::string file_name = "first_mesh_n15";
    std::string file_name = "first_mesh_n50";
    //std::string file_name = "first_mesh_no_hole";

    auto mesh = Geometry::Mesh::parse_mesh("../meshes/" + file_name + ".msh");

    // temporary solution_vector
    std::vector<f32> z;
    for (int i = 0; i < mesh.nof_vertices; ++i) {
        z.push_back(std::abs(mesh.vertices.at(i).at(0) * mesh.vertices.at(i).at(1)));
    }

    mesh.save_mesh_3D(file_name + ".vtk", z);
    
    return 0;
}
