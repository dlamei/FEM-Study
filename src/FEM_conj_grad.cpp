#include "FEM_naive.h"

#include "MatrixImpl/sparse_matrix.h"
#include "MatrixImpl/stack_matrix.h"
#include "MatrixImpl/conj_grad.h"

#include <fstream>
#include <iostream>

using Vector = std::vector<scalar>;

typedef stack_matrix::Vector3 Vec3;
typedef stack_matrix::Matrix3x3 Mat3x3;
typedef stack_matrix::Matrix3x2 Mat3x2;
typedef stack_matrix::Matrix2x3 Triangle;
typedef sparse_matrix::Matrix SparseMatrix;

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
	Mat3x3 grad3x3 = Mat3x3::init(
		1, tri.get(0, 0), tri.get(1, 0),
		1, tri.get(0, 1), tri.get(1, 1),
		1, tri.get(0, 2), tri.get(1, 2)
	);
	grad3x3.inverse_inplace();
	Triangle grad2x3 = Triangle::init( // could be way better if we can convert a 3x3 matrix to a 2x3 matrix by cutting
		grad3x3.get(1, 0), grad3x3.get(1, 1), grad3x3.get(1, 2),
		grad3x3.get(2, 0), grad3x3.get(2, 1), grad3x3.get(2, 2)
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

Mat3x3 local_elem_mat(const Triangle &tri) {
	PROFILE_FUNC()


		scalar area = triangle_area(tri);
	Triangle grad = grad_bary_coords(tri);
	Mat3x2 grad_t = grad.transpose();
	Mat3x3 grad3x3 = stack_matrix::mat3x2_mul_mat2x3(grad_t, grad);
	grad3x3.mul(area);

	return grad3x3;
}


Vec3 local_load_vec(const Triangle &tri) {
	Vec3 local_phi{ 0,0,0 };
	return local_phi;
}


//* ASSEMBLING *//


/*
 * assemble galerkin matrix
 *
 * element / stiffness matrix
 *
 */

std::vector<Triplet> assemble_galerkin_mat(const Mesh &mesh) {
	PROFILE_FUNC()

		usize n_verts = mesh.n_nodes;
	usize n_cells = mesh.n_triangles;

	// TODO: #nnz estimation for resizing triplets vector
	std::vector<Triplet> triplets;

	for (u32 i = 0; i < n_cells; i++) {

		Tria tri_indices = mesh.triangles.at(i);

		Triangle tri = tri_from_indx(mesh, i);
		Mat3x3 elem_mat = local_elem_mat(tri);

		// loop unrolled of
		//
		//for j in  0..3
		//    for k in 0..3

		//if (i == 0) {
		//	std::cout << "Tri Incices: \n";
		//	std::cout << tri_indices.a << " " << tri_indices.b << " " << tri_indices.c << "\n";
		//	std::cout << "Triangle: \n";
		//	tri.print();
		//	std::cout << "\n";
		//	std::cout << "\nGrad Bary Coords: \n";
		//	grad_bary_coords(tri).print();
		//	std::cout << "\nElement Matrix: \n";
		//	elem_mat.print();
		//}

		index_t indx_0 = tri_indices.a;
		index_t indx_1 = tri_indices.b;
		index_t indx_2 = tri_indices.c;

		triplets.push_back({ indx_0, indx_0, elem_mat.get(0, 0) });
		triplets.push_back({ indx_1, indx_0, elem_mat.get(1, 0) });
		triplets.push_back({ indx_2, indx_0, elem_mat.get(2, 0) });
		triplets.push_back({ indx_0, indx_1, elem_mat.get(0, 1) });
		triplets.push_back({ indx_1, indx_1, elem_mat.get(1, 1) });
		triplets.push_back({ indx_2, indx_1, elem_mat.get(2, 1) });
		triplets.push_back({ indx_0, indx_2, elem_mat.get(0, 2) });
		triplets.push_back({ indx_1, indx_2, elem_mat.get(1, 2) });
		triplets.push_back({ indx_2, indx_2, elem_mat.get(2, 2) });
	}

	/* write_sparse_to_file(galerkin, "galerkin_mat.txt"); */

	return triplets;
}



Vector assemble_load_vec(const Mesh &mesh) {
	usize n_nodes = mesh.n_nodes;
	usize n_tris = mesh.n_triangles;

	Vector phi(n_nodes, 0);

	for (usize i = 0; i < n_tris; i++) {
		Tria tri_indices = mesh.triangles.at(i);
		Triangle tri = tri_from_indx(mesh, i);

		Vec3 local_phi = local_load_vec(tri);

		phi.at(tri_indices.a) += local_phi.x;
		phi.at(tri_indices.b) += local_phi.y;
		phi.at(tri_indices.c) += local_phi.z;
	}

	return phi;
}

//* SOLVING *//
//bool isSymmetric(const SparseMatrix &matrix, double &unsymmetryMeasure) {
//	int size = matrix.height;
//	unsymmetryMeasure = 0;
//
//	for (int i = 0; i < size; ++i) {
//		for (int j = i; j < size; ++j) { // Check only half the matrix (above the diagonal)
//			if (matrix.get(i, j) != matrix.get(j, i)) {
//				unsymmetryMeasure += abs(matrix.get(i, j) - matrix.get(j, i));
//			}
//		}
//	}
//
//	return unsymmetryMeasure == 0;
//}

void apply_boundry_conditions(const Mesh &mesh, std::vector<Triplet> *A_triplets, Vector *phi) {

	std::set<int> boundry_nodes;
	std::set<int> inner_boundry_nodes;
	std::set<int> outer_boundry_nodes;
	for (auto row : mesh.inner_boundries) {
		boundry_nodes.insert(row);
		inner_boundry_nodes.insert(row);
	}
	for (auto row : mesh.outer_boundries) {
		boundry_nodes.insert(row);
	}

	// RHS Vector
	for (auto outer_boundry_node : mesh.outer_boundries) {
		phi->at(outer_boundry_node) = 0;
	}
	for (auto inner_boundry_node : mesh.inner_boundries) {
		phi->at(inner_boundry_node) = 1;
	}

	for (auto t : *A_triplets) {
		if (inner_boundry_nodes.find(t.col) != inner_boundry_nodes.end() &&
			boundry_nodes.find(t.row) == boundry_nodes.end()) {
			phi->at(t.row) += t.val;
		}
	}

	// Galerkin Matrix
	A_triplets->erase(std::remove_if(A_triplets->begin(), A_triplets->end(),
		[&boundry_nodes](const Triplet &t) {
			return (boundry_nodes.find(t.row) != boundry_nodes.end()) ||
			(boundry_nodes.find(t.col) != boundry_nodes.end());
		}), A_triplets->end());

	for (auto row : mesh.inner_boundries) { A_triplets->push_back({ row, row, 1.0 }); }
	for (auto row : mesh.outer_boundries) { A_triplets->push_back({ row, row, 1.0 }); }

}

std::vector<scalar> solve_fem(const Mesh &mesh) {

	std::vector<Triplet> A_triplets = assemble_galerkin_mat(mesh);
	Vector phi = assemble_load_vec(mesh);
	apply_boundry_conditions(mesh, &A_triplets, &phi);

	index_t n_verts = mesh.n_nodes;
	SparseMatrix A = SparseMatrix::from_triplets(n_verts, n_verts, &A_triplets);


	return fem_ocl::run_conj_grad(A, phi);
}



// mesh utility functions
inline scalar triangle_area(const Triangle &t) {
	return 0.5 * std::abs((t.get(0, 1) - t.get(0, 0)) * (t.get(1, 2) - t.get(1, 1)) -
		(t.get(0, 2) - t.get(0, 1)) * (t.get(1, 1) - t.get(1, 0)));
}


inline Triangle tri_from_indx(const Mesh &mesh, usize indx) {
	Tria tri_indices = mesh.triangles.at(indx);
	Triangle vert_coords;
	vert_coords.set(0, 0, mesh.vertices.at(tri_indices.a).x);
	vert_coords.set(0, 1, mesh.vertices.at(tri_indices.b).x);
	vert_coords.set(0, 2, mesh.vertices.at(tri_indices.c).x);
	vert_coords.set(1, 0, mesh.vertices.at(tri_indices.a).y);
	vert_coords.set(1, 1, mesh.vertices.at(tri_indices.b).y);
	vert_coords.set(1, 2, mesh.vertices.at(tri_indices.c).y);

	return vert_coords;
}


//void write_sparse_to_file(const SparseMatrix &m, const char *name) {
//	PROFILE_FUNC()
//
//		std::ofstream file(name);
//	if (file.is_open())
//	{
//		for (int i = 0; i < m.height; i++) {
//			for (int j = 0; j < m.width; j++) {
//				file << m.get(i, j);
//				file << " ";
//			}
//			file << "\n";
//		}
//	}
//}
