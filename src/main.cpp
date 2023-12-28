#include <iostream>
#include "utils.h"
#include "./Geometries/geometry.h"
#include <cmath>

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
    //std::string file_name = "first_mesh_n50";
    //std::string file_name = "first_mesh_n100";
    //std::string file_name = "first_mesh_n500";
    //std::string file_name = "first_mesh_no_hole";
    std::string file_name = "4_holes_mesh";

    auto mesh = Geometry::Mesh::parse_mesh("../meshes/" + file_name + ".msh");

    // temporary solution_vector
    std::vector<Geometry::scalar_t> z;
    for (int i = 0; i < mesh.nof_vertices; ++i) {
        Geometry::scalar_t x = mesh.vertices.at(i).at(0);
        Geometry::scalar_t y = mesh.vertices.at(i).at(1);
        z.push_back((std::exp(-(x - 0.75) * (x - 0.75)) * std::exp(-(y + 0.5) * (y + 0.5))));
    }

    mesh.save_mesh_3D(file_name + ".vtk", z);
    
    return 0;
}
