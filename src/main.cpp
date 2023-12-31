#include "matrix.h"

//template <typename T> void print_mat(std::vector<std::vector<T>> mat) {
//  for (int i = 0; i < mat.size(); ++i) {
//    for (int j = 0; j < mat.at(0).size(); ++j) {
//      std::cout << mat.at(i).at(j) << " ";
//    }
//    std::cout << "\n";
//  }
//}

int main() {

	auto m1 = matrix::ident(3, 3);
	auto m2 = matrix::ident(3, 3);
	matrix::print(matrix::add(m1, m2));

	/* auto m2 = matrix::ident(2, 2); */

	/* matrix::print(m1); */

  // std::string file_name = "first_mesh_n15";
  // std::string file_name = "first_mesh_n50";
  // std::string file_name = "first_mesh_n100";
  // std::string file_name = "first_mesh_n500";
  // std::string file_name = "first_mesh_no_hole";
  //std::string file_name = "4_holes_mesh";

  //auto mesh = Geometry::Mesh::parse_mesh("../meshes/" + file_name + ".msh");

  //std::vector<scalar> z;
  //for (int i = 0; i < mesh.nof_vertices; ++i) {
  //  scalar x = mesh.vertices.at(i).at(0);
  //  scalar y = mesh.vertices.at(i).at(1);
  //  z.push_back((std::exp(-(x - 0.75) * (x - 0.75)) *
  //               std::exp(-(y + 0.5) * (y + 0.5))));
  //}

  //mesh.save_mesh_3D(file_name + ".vtk", z);

  return 0;
}
