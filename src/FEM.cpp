#include "FEM.h"
#include <fstream>
#include <iostream>

const int Dynamic = Eigen::Dynamic;

template <const int rows, const int cols>
using Matrix = Eigen::Matrix<scalar, rows, cols>;

template <const int rows>
using Vector = Eigen::Matrix<scalar, rows, 1>;

typedef Matrix<2, 3> Triangle;
typedef Eigen::SparseMatrix<scalar> SparseMatrix;

struct Mesh {

    usize n_nodes { 0 };
    usize n_triangles { 0 };
    usize n_boundary_edges { 0 };

    //TODO: Samuel add boundary cond.

    // 2 scalar per node
    Matrix<Dynamic, 2> nodes {};

    // 3 vert_indx's per triangle
    Eigen::Matrix<index_t , Eigen::Dynamic, 3> triangles{};

};



scalar triangle_area(const Triangle &t);
Triangle tri_from_indx(const Mesh &mesh, usize indx);

void write_sparse_to_file(const SparseMatrix &m, const char *name);

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

Triangle grad_bary_coords(const Triangle &tri) {
    Matrix<3, 3> grad;
    grad.block<3, 1>(0, 0) = Matrix<3, 1>::Ones();
    grad.block<3, 2>(0, 1) = tri.transpose();
    return grad.inverse().block<2, 3>(1, 0);
}


//* LOCAL COMPUTATIONS *//

/*
 * local galerkin matrix of triangle K
 *
 * a_K(b(j), b(i))
 * 
 * a_K = |K| grad^T * grad
 *
 */

Matrix<3, 3> local_elem_mat(const Triangle &tri) {
    PROFILE_FUNC()

    scalar area = triangle_area(tri);
    auto grad = grad_bary_coords(tri);
    return area * grad.transpose() * grad;
}


inline Vector<3> local_load_vec(const Triangle &tri, source_fn_ptr source_fn) {
    scalar area = triangle_area(tri);

    Vector<3> local_phi;
    local_phi(0) = source_fn(tri.col(0));
    local_phi(1) = source_fn(tri.col(1));
    local_phi(2) = source_fn(tri.col(2));

    local_phi *= area / 3.0;
    return local_phi;
}


//* ASSEMBLING *//


/*
 * assemble galerkin matrix 
 *
 * element / stiffness matrix
 *
 */

SparseMatrix assemble_galerkin_mat(const Mesh &mesh) {
    PROFILE_FUNC()

    usize n_verts = mesh.n_nodes;
    usize n_cells = mesh.n_triangles;

    // TODO: #nnz estimation for resizing triplets vector
    std::vector<Eigen::Triplet<scalar>> triplets;

    for (u32 i = 0; i < n_cells; i++) {

        Eigen::Vector<index_t, 3> tri_indices = mesh.triangles.row(i);

        Triangle tri = tri_from_indx(mesh, i);
        Matrix<3, 3> elem_mat = local_elem_mat(tri);

        // loop unrolled of
        //
        //for j in  0..3
        //    for k in 0..3

        index_t indx_0 = tri_indices(0);
        index_t indx_1 = tri_indices(1);
        index_t indx_2 = tri_indices(2);

        triplets.push_back({ indx_0, indx_0, elem_mat(0, 0) });
        triplets.push_back({ indx_1, indx_0, elem_mat(1, 0) });
        triplets.push_back({ indx_2, indx_0, elem_mat(2, 0) });

        triplets.push_back({ indx_0, indx_1, elem_mat(0, 1) });
        triplets.push_back({ indx_1, indx_1, elem_mat(1, 1) });
        triplets.push_back({ indx_2, indx_1, elem_mat(2, 1) });

        triplets.push_back({ indx_0, indx_2, elem_mat(0, 2) });
        triplets.push_back({ indx_1, indx_2, elem_mat(1, 2) });
        triplets.push_back({ indx_2, indx_2, elem_mat(2, 2) });
    }

    SparseMatrix galerkin(n_verts, n_verts);
    galerkin.setFromTriplets(triplets.begin(), triplets.end());

    /* write_sparse_to_file(galerkin, "galerkin_mat.txt"); */

    return galerkin;
}



Vector<Dynamic> assemble_load_vec(const Mesh &mesh, source_fn_ptr source_fn) {
    usize n_nodes = mesh.n_nodes;
    usize n_tris = mesh.n_triangles;

    Vector<Dynamic> phi = Vector<Dynamic>::Zero(n_nodes);

    for (usize i = 0; i < n_tris; i++) {
        Eigen::Vector<index_t, 3> tri_indices = mesh.triangles.row(i);
        Triangle tri = tri_from_indx(mesh, i);

        Vector<3> local_phi = local_load_vec(tri, source_fn);

        phi(tri_indices(0)) += local_phi(0);
        phi(tri_indices(1)) += local_phi(1);
        phi(tri_indices(1)) += local_phi(2);
    }

    return phi;
}

Mesh geometry_to_mesh(const Geometry &geo) {

    auto nodes = Matrix<Dynamic, 2>::Zero(geo.nof_vertices, 2).eval();
    auto triangles = Eigen::Matrix<index_t, Dynamic, 3>::Zero(geo.nof_triangles, 3).eval();

    for (index_t i = 0; i < geo.vertices.size(); i++) {
        const auto &vert = geo.vertices.at(i);
        nodes(i, 0) = vert.at(0);
        nodes(i, 1) = vert.at(1);
    }

    for (index_t i = 0; i < geo.triangles.size(); i++) {
        const auto &tri = geo.triangles.at(i);
        triangles(i, 0) = tri.at(0);
        triangles(i, 1) = tri.at(1);
        triangles(i, 2) = tri.at(2);
    }

    return Mesh {
        .n_nodes = geo.nof_vertices,
        .n_triangles = geo.nof_triangles,
        .n_boundary_edges = geo.nof_boundry_edges,
        .nodes = nodes,
        .triangles = triangles,
    };
}

//* SOLVING *//

Vector<Dynamic> solve_fem(const Geometry &geo, source_fn_ptr source_fn) {
    auto mesh = geometry_to_mesh(geo);

    SparseMatrix a = assemble_galerkin_mat(mesh);
    Vector<Dynamic> phi = assemble_load_vec(mesh, source_fn);

    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(a);
    solver.factorize(a);
    Vector<Dynamic> mu = solver.solve(phi);

    return mu;
}



// mesh utility functions

inline scalar triangle_area(const Triangle &t) {
    return 0.5 * std::abs((t(0, 1) - t(0, 0)) * (t(1, 2) - t(1, 1)) -
            (t(0, 2) - t(0, 1)) * (t(1, 1) - t(1, 0)));
}


inline Triangle tri_from_indx(const Mesh &mesh, usize indx) {
    Eigen::Vector<index_t, 3> tri_indices = mesh.triangles.row(indx);

    Matrix<3, 2> vert_coords;
    vert_coords << mesh.nodes.row(tri_indices[0]),
                mesh.nodes.row(tri_indices[1]),
                mesh.nodes.row(tri_indices[2]);

    return vert_coords.transpose();
}


void write_sparse_to_file(const SparseMatrix &m, const char *name) {
    PROFILE_FUNC()

    std::ofstream file(name);
    if (file.is_open())
    {
        file << m << std::endl;
    }
}
