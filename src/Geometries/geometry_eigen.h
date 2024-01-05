#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include "../utils.h"

template <const int rows, const int cols>
using Matrix = Eigen::Matrix<scalar, rows, cols>;

template <const int rows>
using Vector = Eigen::Matrix<scalar, rows, 1>;

typedef Matrix<2, 3> Triangle;
typedef Eigen::SparseMatrix<scalar> SparseMatrix;

typedef scalar (*source_fn_ptr)(const Vector<2> &);

namespace Geometry {

struct Mesh {
    usize nof_nodes{ 0 };
    usize nof_triangles{ 0 };
    usize nof_boundry_edges{ 0 };

    // Matrix containing x and y coordinates of nodes
    // dimension: nof_nodes * 2
    Matrix<Eigen::Dynamic, 2> nodes{};

    // Matrix containing indecies of the nodes of the triangles
    // dimension: nof_triangles * 3
    Eigen::Matrix<int, Eigen::Dynamic, 3> triangles{};
    
    // Vector containing the nodes which are part of the different boundries
    // boundries.at(0) contains the verticies of all boundries
    // dimension: (nof_different_boundries + 1) * nof_specific boundries * 2
    std::vector<std::vector<std::vector<u64>>> boundries;

    // parses a .msh file into a Mesh object and returns it
    static Mesh parse_mesh(std::string file_name);

    // returns a 2*3 matrix that contains the coordinates of the nodes of the triangele with index i
    Triangle get_tria_coords(int i) const;

    // creates a .vtk file of the solution 
    void save_mesh_3D(std::string filename, const Vector<Eigen::Dynamic> &z) const;

};

}
