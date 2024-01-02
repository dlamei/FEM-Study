#pragma once

#include "../utils.h"
#include <vector>
#include <string>
#include <fstream>

namespace Geometry {
    
using matrix_t = std::vector<std::vector<scalar>>;
using vector_t = std::vector<scalar>;

struct Mesh {
    usize nof_vertices{};
    usize nof_triangles{};
    usize nof_boundry_edges{};

    // Matrix containing x and y coordinates of vertices
    // dimension: nof_vertices * 2
    matrix_t vertices;

    // Matrix containing indecies of the vertices of the triangles
    // dimension: nof_triangles * 3
    std::vector<std::vector<u64>> triangles;
    
    // Vector containing the vertices which are part of the different boundries
    // boundries.at(0) contains the verticies of all boundries
    // dimension: (nof_different_boundries + 1) * nof_specific boundries * 2
    std::vector<std::vector<std::vector<u64>>> boundries;

    // parses a .msh file into a Mesh object and returns it
    static Mesh parse_mesh(std::string file_name);

    // returns a 2*3 matrix that contains the coordinates of the vertices of the triangele with index i
    matrix_t get_tria_coords(int i);

    // creates a .vtk file of the solution 
    void save_mesh_3D(std::string filename, const std::vector<scalar> &z) const;

};

}
