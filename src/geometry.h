#pragma once

#include "utils.h"
#include <vector>
#include <string>
#include <fstream>

//* generic way to parse and store mesh files *//
struct Coords {
    scalar x{ 0 };
    scalar y{ 0 } ;
};
struct Tria {
    index_t a{ 0 };
    index_t b{ 0 };
    index_t c{ 0 };
};

struct Mesh {
    usize n_nodes{ 0 };
    usize n_triangles{ 0 };
    usize n_boundry_nodes{ 0 };

    // Matrix containing x and y coordinates of vertices
    // dimension: nof_vertices * 2
    std::vector<Coords> vertices;

    // Matrix containing indecies of the vertices of the triangles
    // dimension: nof_triangles * 3
    std::vector<Tria> triangles;
    
    // Vector containing the vertices which are part of the inner boundry
    // dimension: nof_inner_boundry_nodes * 1
    std::vector<index_t> inner_boundries;

    // Vector containing the vertices which are part of the outer boundry
    // dimension: nof_outer_boundry_nodes * 1
    std::vector<index_t> outer_boundries;

    // parses a .msh file into a Mesh object and returns it
    static Mesh parse_mesh(std::string file_name);

    // creates a .vtk file of the solution
    // takes as input a pointer to the resulting vector, count is the # of elements in result
    // every implementation returns different type of vector as result
    void save_mesh_3D(const std::string& filename, const scalar *result, usize count) const;
};
