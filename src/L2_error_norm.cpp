#include <chrono>
#include <fstream>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include "FEM_naive.h"

#define MAX_N 400

int main() {

    std::chrono::time_point< std::chrono::high_resolution_clock > t_start, t_end;

    std::ofstream output("time_series.txt");
    assert(output.is_open(),
            "Failed to open the mesh file");
    
    std::string file_name = "n_1h";
    std::string folder_name = "1h_linear_series/";

    Eigen::Vector2d a(-16.0 / 15.0, 0.0);
    Eigen::Vector2d b(-1.0 / 15.0, 0.0);
    // Lambda function providing exact solution
    auto uExact = [&a, &b](Eigen::VectorXd x) -> double {
      return (log((x - a).norm()) - log((x - b).norm())) / log(2) - 1;
    };

    for(int i = 10; i <= MAX_N; i = i + 10) {
        auto mesh = Mesh::parse_mesh("../meshes/" + folder_name + std::to_string(i) + file_name + ".msh");
        auto res = solve_fem(mesh);
        double errorSum = 0.0;
        for (int j = 0; j < mesh.n_nodes; ++j) {
            Eigen::Vector2d point(); // Construct the point (x, 0)
            double exactValue = uExact(point); // Evaluate the exact solution
            double error = exactValue - approx[i]; // Compute the difference
            errorSum += error * error; // Add the square of the difference to the sum
        }

    }


    return 0;
}
