#include "FEM_naive.h"
#include <fstream>
#include <iostream>


using Matrix = naive_matrix::Matrix;

using Vector = naive_matrix::Matrix;

typedef stack_matrix::Vector3 Vector3;
typedef stack_matrix::Matrix3x3 Matrix3x3;
typedef stack_matrix::Matrix2x3 Triangle;
typedef naive_matrix::Matrix SparseMatrix;

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
    auto grad3x3 = stack_matrix::Matrix3x3::init(
                                    1, tri.get(0,0), tri.get(0,1),
                                    1, tri.get(1,0), tri.get(1,1),
                                    1, tri.get(2,0), tri.get(2,1)
                                    );
    grad3x3.inverse_inplace();
    auto grad2x3 = stack_matrix::Matrix2x3::init( // could be way better if we can convert a 3x3 matrix to a 2x3 matrix by cutting
                                    grad3x3.get(1,0), grad3x3.get(1,1), grad3x3.get(1,2),
                                    grad3x3.get(2,0), grad3x3.get(2,1), grad3x3.get(2,2)
                                    );
    return grad2x3;
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

Matrix3x3 local_elem_mat(const Triangle &tri) {
    PROFILE_FUNC()

    scalar area = triangle_area(tri);
    auto grad = grad_bary_coords(tri);
    auto grad_t = grad.transpose();
    auto grad3x3 = stack_matrix::mat3x2_mul_mat2x3(grad_t, grad);
    grad3x3.mul(area);
    
    return grad3x3;
}


Vector3 local_load_vec(const Triangle &tri) {
    Vector3 local_phi{0,0,0};
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
    SparseMatrix gal_matrix = naive_matrix::Matrix::zero(n_cells, n_cells);

    for (u32 i = 0; i < n_cells; i++) {

        auto tri_indices = mesh.triangles.at(i);

        Triangle tri = tri_from_indx(mesh, i);
        Matrix3x3 elem_mat = local_elem_mat(tri);

        // loop unrolled of
        //
        //for j in  0..3
        //    for k in 0..3

        index_t indx_0 = tri_indices.a;
        index_t indx_1 = tri_indices.b;
        index_t indx_2 = tri_indices.c;

        gal_matrix.set( indx_0, indx_0, elem_mat.get(0, 0) );
        gal_matrix.set( indx_1, indx_0, elem_mat.get(1, 0) );
        gal_matrix.set( indx_2, indx_0, elem_mat.get(2, 0) );
        gal_matrix.set( indx_0, indx_1, elem_mat.get(0, 1) );
        gal_matrix.set( indx_1, indx_1, elem_mat.get(1, 1) );
        gal_matrix.set( indx_2, indx_1, elem_mat.get(2, 1) );
        gal_matrix.set( indx_0, indx_2, elem_mat.get(0, 2) );
        gal_matrix.set( indx_1, indx_2, elem_mat.get(1, 2) );
        gal_matrix.set( indx_2, indx_2, elem_mat.get(2, 2) );
    }

    // Boundry conditions
    std::set<int> rowsToRemove;
    for(auto row : mesh.inner_boundries) {
        rowsToRemove.insert(row);
    }
    for(auto row : mesh.outer_boundries) {
        rowsToRemove.insert(row);
    }

    triplets.erase(std::remove_if(triplets.begin(), triplets.end(), 
        [&rowsToRemove](const Eigen::Triplet<scalar> &t) {
            return rowsToRemove.find(t.row()) != rowsToRemove.end();
        }), triplets.end());

    for(auto row : mesh.inner_boundries) {
        triplets.push_back(Eigen::Triplet<scalar>(row, row, 1.0));
    }
    for(auto row : mesh.outer_boundries) {
        triplets.push_back(Eigen::Triplet<scalar>(row, row, 1.0));
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
        auto tri_indices = mesh.triangles.at(i);
        Triangle tri = tri_from_indx(mesh, i);

        Vector<3> local_phi = local_load_vec(tri, source_fn);

        phi(tri_indices.a) += local_phi(0);
        phi(tri_indices.b) += local_phi(1);
        phi(tri_indices.c) += local_phi(2);
    }

    // Boundry Conditions
    for(auto outer_boundry_node : mesh.outer_boundries) {
        phi(outer_boundry_node) = 0.0   ;
    }
    for(auto inner_boundry_node : mesh.inner_boundries) {
        phi(inner_boundry_node) = 1.0;
    }

    return phi;
}

//* SOLVING *//

naive_matrix::Matrix solve_fem(const Mesh &mesh) {
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
    return 0.5 * std::abs((t.get(0, 1) - t.get(0, 0)) * (t.get(1, 2) - t.get(1, 1)) -
            (t.get(0, 2) - t.get(0, 1)) * (t.get(1, 1) - t.get(1, 0)));
}


inline Triangle tri_from_indx(const Mesh &mesh, usize indx) {
    auto tri_indices = mesh.triangles.at(indx);
    Matrix<3, 2> vert_coords;
    vert_coords << mesh.vertices.at(tri_indices.a).x,
                   mesh.vertices.at(tri_indices.b).x,
                   mesh.vertices.at(tri_indices.c).x,
                   mesh.vertices.at(tri_indices.a).y,
                   mesh.vertices.at(tri_indices.b).y,
                   mesh.vertices.at(tri_indices.c).y;

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
