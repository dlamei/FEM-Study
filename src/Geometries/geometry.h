#pragma once

#include "../utils.h"
#include <vector>
#include <string>
#include <fstream>

namespace Geometry {

using scalar_t = f64;

struct Mesh {
    usize nof_vertices{};
    usize nof_triangles{};
    usize nof_boundry_edges{};

    // Matrix containing x and y coordinates of vertices
    // dimension: nof_vertices * 2
    std::vector<std::vector<scalar_t>> vertices;

    // Matrix containing indecies of the vertices of the triangles
    // dimension: nof_triangles * 3
    std::vector<std::vector<u64>> triangles;
    
    // Vector containing the vertices wich are on the outer boundry
    std::vector<std::vector<u64>> outer_boundry;

    // Vector containing the vertices wich are on the  inner boundry
    std::vector<std::vector<u64>> inner_boundry;

    // parses a .msh file into a Mesh object and returns it
    static Mesh parse_mesh(std::string file_name) {
        Mesh mesh{};

        std::fstream input;
        input.open(file_name);
        assert(input.is_open(), 
                "Failed to open the mesh file");
        
        // read in number of vertices, triangeles and edge vertices
        input >> mesh.nof_vertices;
        input >> mesh.nof_triangles;
        input >> mesh.nof_boundry_edges;

        assert((mesh.nof_vertices >= 0 && mesh.nof_triangles >= 0 && mesh.nof_boundry_edges >= 0), 
                "Error nof_vertices, nof_triangles or nof_boundry_edges is an invalid input");
        
        // read in vertice coordinates
        for(int i = 0; i < mesh.nof_vertices; ++i) {
            scalar_t x, y;
            u64 boundry_label;
            input >> x;
            input >> y;
            input >> boundry_label;
            mesh.vertices.push_back({x,y});
        }

        // read in the triangles
        for(int i = 0; i < mesh.nof_triangles; ++i) {
            u64 l, m, n;
            u64 boundry_label;
            input >> l;
            input >> m;
            input >> n;
            input >> boundry_label;
            mesh.triangles.push_back({l, m, n});
        }

        // read in the boundries
        for(int i = 0; i < mesh.nof_boundry_edges; ++i) {
            u64 m, n;
            u64 boundry_label;
            input >> m;
            input >> n;
            input >> boundry_label;
            if(boundry_label == 1) {
                mesh.outer_boundry.push_back({m,n});
            }
            else if(boundry_label == 2) {
                mesh.inner_boundry.push_back({m,n});
            }
        }

        return mesh;
    }

    // returns a 2*3 matrix that contains the coordinates of the vertices of the triangele with index i
    std::vector<std::vector<scalar_t>> get_tria_coords(int i) {
        std::vector<std::vector<scalar_t>> tria_coords;
        tria_coords.push_back(vertices.at(triangles.at(i).at(0) - 1));
        tria_coords.push_back(vertices.at(triangles.at(i).at(1) - 1));
        tria_coords.push_back(vertices.at(triangles.at(i).at(2) - 1));
        return tria_coords;
    }

};

}