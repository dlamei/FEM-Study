#include "FEM.h"
#include <fstream>
#include <iostream>




scalar triangle_area(const Triangle &t);
Triangle tri_from_indx(const Mesh &mesh, usize indx);

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

SparseMatrix assemble_gelerkin_mat(const Mesh &mesh) {

    usize n_verts = mesh.n_nodes;
    usize n_cells = mesh.n_triangles;

    // TODO: #nnz estimation for resizing triplets vector
    std::vector<Eigen::Triplet<scalar>> triplets;

    for (u32 i = 0; i < n_cells; i++) {

        Eigen::Vector3i tri_indices = mesh.triangles.row(i);

        Triangle tri = tri_from_indx(mesh, i);
        Matrix<3, 3> elem_mat = local_elem_mat(tri);

        // loop unrolled of
        //
        //for j in  0..3
        //    for k in 0..3

        i32 indx_0 = tri_indices(0);
        i32 indx_1 = tri_indices(1);
        i32 indx_2 = tri_indices(2);

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

    return galerkin;
}



Vector<Dim::Dynamic> assemble_load_vec(const Mesh &mesh, source_fn_ptr source_fn) {
    usize n_nodes = mesh.n_nodes;
    usize n_tris = mesh.n_triangles;

    Vector<Dim::Dynamic> phi = Vector<Dim::Dynamic>::Zero(n_nodes);

    for (usize i = 0; i < n_tris; i++) {
        Eigen::Vector3i tri_indices = mesh.triangles.row(i);
        Triangle tri = tri_from_indx(mesh, i);

        Vector<3> local_phi = local_load_vec(tri, source_fn);

        phi(tri_indices(0)) += local_phi(0);
        phi(tri_indices(1)) += local_phi(1);
        phi(tri_indices(1)) += local_phi(2);
    }

    return phi;
}

//* SOLVING *//

Vector<Dim::Dynamic> solve_fem(const Mesh &mesh, source_fn_ptr source_fn) {
    SparseMatrix a = assemble_gelerkin_mat(mesh);
    Vector<Dim::Dynamic> phi = assemble_load_vec(mesh, source_fn);

    Eigen::SparseLU<SparseMatrix, Eigen::COLAMDOrdering<int>> solver;
    solver.analyzePattern(a);
    solver.factorize(a);
    Vector<Dim::Dynamic> mu = solver.solve(phi);

    return mu;
}



// mesh utility functions

inline scalar triangle_area(const Triangle &t) {
    return 0.5 * std::abs((t(0, 1) - t(0, 0)) * (t(1, 2) - t(1, 1)) -
            (t(0, 2) - t(0, 1)) * (t(1, 1) - t(1, 0)));
}


inline Triangle tri_from_indx(const Mesh &mesh, usize indx) {
    Eigen::Vector3i tri_indices = mesh.triangles.row(indx);

    Matrix<3, 2> vert_coords;
    vert_coords << mesh.nodes.row(tri_indices[0]),
                mesh.nodes.row(tri_indices[1]),
                mesh.nodes.row(tri_indices[2]);

    return vert_coords.transpose();
}


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

    auto nodes = Eigen::Matrix<scalar, Dim::Dynamic, 2>(n_verts, 2);
    for (usize i = 0; i < n_verts; i++) {
        scalar x{}, y{};
        u64 boundary_label;
        input >> x;
        input >> y;
        input >> boundary_label;
        nodes(i, 0) = x;
        nodes(i, 1) = y;
    }

    auto tris = Eigen::Matrix<int, Dim::Dynamic, 3>(n_tris, 3);
    for(int i = 0; i < n_tris; ++i) {
        u64 l, m, n;
        u64 boundry_label;
        input >> l;
        input >> m;
        input >> n;
        input >> boundry_label;
        tris(i, 0) = l - 1;
        tris(i, 1) = m - 1;
        tris(i, 2) = n - 1;
    }


    Mesh mesh {
        .n_nodes = n_verts,
            .n_triangles = n_tris,
            .nodes = nodes,
            .triangles = tris,
    };

    return mesh;

};
