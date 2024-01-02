#include "FEM.h"
#include <fstream>
#include <iostream>

Mesh Mesh::load(const std::string &file_name) {
    std::fstream input;
    input.open(file_name);

    assert(input.is_open(), "Failed to open mesh file");

    usize n_verts = 0;
    usize n_tris = 0;
    usize n_boundry_edges = 0;
    input >> n_verts;
    input >> n_tris;
    input >> n_boundry_edges;
    assert((n_verts >= 0 && n_tris >= 0 && n_boundry_edges >= 0), 
            "Error nof_vertices, nof_triangles or nof_boundry_edges is an invalid input");

    auto nodes = Eigen::Matrix<scalar, Eigen::Dynamic, 2>(n_verts, 2);
    for (usize i = 0; i < n_verts; i++) {
        scalar x{}, y{};
        u64 boundary_label;
        input >> x;
        input >> y;
        input >> boundary_label;
        nodes(i, 0) = x;
        nodes(i, 1) = y;
    }

    auto elems = Eigen::Matrix<int, Eigen::Dynamic, 3>(n_tris, 3);
    for(int i = 0; i < n_tris; ++i) {
        u64 l, m, n;
        u64 boundry_label;
        input >> l;
        input >> m;
        input >> n;
        input >> boundry_label;
        elems(i, 0) = l;
        elems(i, 1) = m;
        elems(i, 2) = n;
    }


    Mesh mesh {
        .nof_vertices = n_verts,
            .nof_triangles = n_tris,
            .nof_boundary_edges = n_boundry_edges,
            .nodes = nodes,
            .triangles = elems,
    };

    return mesh;

};

Mesh::TriGeo Mesh::get_tria_coords(usize i) {
    assert(i  < triangles.rows());

    const auto indx = triangles.row(i);
    Eigen::Matrix<scalar, 3, 2> verts;
    verts << nodes.row(indx[0]),
          nodes.row(indx[1]),
          nodes.row(indx[2]);

    return verts.transpose();
};


/*
 *  l1, l2, l3 on triangle K => barycentric coords functions
 *                                                           
 *  1 a(1, 1) a(1, 2)   y(1)    y(2)    y(3)      1 0 0
 *  1 a(2, 1) a(2, 2) * x(1, 1) x(2, 1) x(3, 1) = 0 1 0
 *  1 a(3, 1) a(3, 2)   x(1, 2) x(2, 2) x(3, 2)   0 0 1
 *                                                           
 *  grad (l(i)) = x(i)
 *                                                           
 *  solve for x
 *
 */
 
Mesh::TriGeo grad_bary_coords(const Mesh::TriGeo &verts) {
    Eigen::Matrix<scalar, 3, 3> grad;
    grad.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
    grad.block<3, 2>(0, 1) = verts.transpose();
    return grad.inverse().block<2, 3>(1, 0);
}

/*
 * local galerkin matrix element
 *
 * a_K(b(j), b(i))
 *
 * a_K(1), a_K(2), a_K(3) verts of triangle K
 *
 * 
 * a_K = |K| grad^T * grad
 *
 */

Eigen::Matrix<scalar, 3, 3> elem_matrix(const Mesh::TriGeo &v) {
    scalar area = 0.5 * std::abs((v(0, 1) - v(0, 0)) * (v(1, 2) - v(1, 1)) -
                                 (v(0, 2) - v(0, 1)) * (v(1, 1) - v(1, 0)));

    auto grad = grad_bary_coords(v);

    return area * grad.transpose() * grad;
}

//TODO assemble galerkin + rhs (cell oriented) 2.4.5.24
