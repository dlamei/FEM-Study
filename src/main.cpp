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


    defer(print(0));
    defer(print(1));
    defer(print(2));
    defer(print(3));

    std::string file_name = "first_mesh_n15.msh";
    auto geo = Geometry::Mesh::parse_mesh("../meshes/" + file_name);

    for (int i = 0; i < 5; ++i) {
        print_mat(geo.get_tria_coords(0));
        std::cout << "\n";
    }
    
    assert(1 == 2, "this is false");

    return 0;
}
