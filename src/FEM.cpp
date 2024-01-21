#include "FEM.h"
#include <fstream>
#include <iostream>


#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>


const int Dynamic = Eigen::Dynamic;

template <const int rows, const int cols>
using Matrix = Eigen::Matrix<scalar, rows, cols>;

template <const int rows>
using Vector = Eigen::Matrix<scalar, rows, 1>;

typedef Matrix<2, 3> Triangle;
typedef Eigen::SparseMatrix<scalar> SparseMatrix;

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
	PROFILE_FUNC();

	usize n_verts = mesh.n_nodes;
	usize n_cells = mesh.n_triangles;

	// TODO: #nnz estimation for resizing triplets vector
	std::vector<Eigen::Triplet<scalar>> triplets;

	for (u32 i = 0; i < n_cells; i++) {

		auto tri_indices = mesh.triangles.at(i);

		Triangle tri = tri_from_indx(mesh, i);
		Matrix<3, 3> elem_mat = local_elem_mat(tri);

		// loop unrolled of
		//
		//for j in  0..3
		//    for k in 0..3

		index_t indx_0 = tri_indices.a;
		index_t indx_1 = tri_indices.b;
		index_t indx_2 = tri_indices.c;

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



	// Boundry conditions
	std::set<int> rowsToRemove;
	for (auto row : mesh.inner_boundries) {
		rowsToRemove.insert(row);
	}
	for (auto row : mesh.outer_boundries) {
		rowsToRemove.insert(row);
	}

	triplets.erase(std::remove_if(triplets.begin(), triplets.end(),
		[&rowsToRemove](const Eigen::Triplet<scalar> &t) {
			return rowsToRemove.find(t.row()) != rowsToRemove.end();
		}), triplets.end());

	for (auto row : mesh.inner_boundries) {
		triplets.push_back(Eigen::Triplet<scalar>(row, row, 1.0));
	}
	for (auto row : mesh.outer_boundries) {
		triplets.push_back(Eigen::Triplet<scalar>(row, row, 1.0));
	}

	SparseMatrix galerkin(n_verts, n_verts);
	galerkin.setFromTriplets(triplets.begin(), triplets.end());
	galerkin.makeCompressed();

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
	for (auto outer_boundry_node : mesh.outer_boundries) {
		phi(outer_boundry_node) = 0.0;
	}
	for (auto inner_boundry_node : mesh.inner_boundries) {
		phi(inner_boundry_node) = 1.0;
	}

	return phi;
}

//void calc_band_width(const SparseMatrix &_m) {
//	auto m = _m.toDense();
//
//	i32 max_width = 0;
//
//	for (i32 i = 0; i < m.rows(); i++) {
//		i32 start = 0;
//		i32 end = 0;
//
//		for (i32 j = 0; j < m.cols(); j++) {
//			scalar v = m(i, j);
//
//			if (!cmp_scalar(v, 0)) {
//				if (start == 0) start = j;
//				end = j;
//			}
//
//		}
//
//		i32 width = end - start;
//		max_width = std::max(max_width, width);
//	}
//
//	printf("band_width: %i\n", max_width);
//}

//* SOLVING *//

Vector<Dynamic> solve_fem(const Mesh &mesh, source_fn_ptr source_fn) {
	SparseMatrix a = assemble_galerkin_mat(mesh);
	Vector<Dynamic> phi = assemble_load_vec(mesh, source_fn);

	Eigen::SparseLU<SparseMatrix, Eigen::NaturalOrdering<i32>> solver;

	solver.analyzePattern(a);
	solver.factorize(a);

	if (solver.info() != Eigen::ComputationInfo::Success) {
		printf("solve: ComputationInfo %d", solver.info());
		exit(-1);
	}

	Vector<Dynamic> mu = solver.solve(phi);


	/* write_sparse_to_file(a, "m.txt"); */
	/* write_sparse_to_file(solver.matrixL().toSparse(), "l.txt"); */
	/* write_sparse_to_file(solver.matrixU().toSparse(), "u.txt"); */



	return mu;
}



// mesh utility functions

inline scalar triangle_area(const Triangle &t) {
	return 0.5 * std::abs((t(0, 1) - t(0, 0)) * (t(1, 2) - t(1, 1)) -
		(t(0, 2) - t(0, 1)) * (t(1, 1) - t(1, 0)));
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
	std::ofstream file(name);

	file << m.rows() << " " << m.cols() << "\n";
	for (int k = 0; k < m.outerSize(); ++k) {
		for (SparseMatrix::InnerIterator it(m, k); it; ++it) {
			file << it.col() << " " << it.row() << " " << it.value() << "\n";
		}
	}
}
