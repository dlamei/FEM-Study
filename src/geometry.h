#pragma once

#include "utils.h"
#include <vector>
#include <string>
#include <fstream>

//* generic way to parse and store mesh files *//

struct Geometry {
    usize nof_vertices{ 0 };
    usize nof_triangles{ 0 };
    usize nof_boundry_edges{ 0 };

    // Matrix containing x and y coordinates of vertices
    // dimension: nof_vertices * 2
    std::vector<std::vector<scalar>> vertices;

    // Matrix containing indecies of the vertices of the triangles
    // dimension: nof_triangles * 3
    std::vector<std::vector<index_t>> triangles;
    
    // Vector containing the vertices which are part of the different boundries
    // boundries.at(0) contains the verticies of all boundries
    // dimension: (nof_different_boundries + 1) * nof_specific boundries * 2
    std::vector<std::vector<std::vector<index_t>>> boundries;

    // parses a .msh file into a Mesh object and returns it
    static Geometry parse_mesh(std::string file_name);

    // creates a .vtk file of the solution
    // takes as input a pointer to the resulting vector, count is the # of elements in result
    // every implementation returns different type of vector as result
    void save_mesh_3D(const std::string& filename, const scalar *result, usize count) const;
};
