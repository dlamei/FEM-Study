#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include "utils.h"

typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;

struct Mesh {

    typedef Eigen::Matrix<scalar, 2, 3> TriGeo;

    usize nof_vertices { 0 };
    usize nof_triangles { 0 };
    usize nof_boundary_edges { 0 };

    Eigen::Matrix<scalar, Eigen::Dynamic, 2> nodes{};
    Eigen::Matrix<int, Eigen::Dynamic, 3> triangles{};

    // get coords of verts of triangle i
    TriGeo get_tria_coords(usize i);

    static Mesh load(const std::string &file_name);
};

