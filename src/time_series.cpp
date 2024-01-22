#include <chrono>
#include <fstream>
#include "FEM_naive.h"



#define MAX_N 80

/*
scalar source_fn(const Eigen::Vector<scalar, 2> &x) {
    return 0;
}
*/


// Function to compute the Euclidean norm of the difference of two 2D vectors
scalar norm(const std::array<scalar, 2>& vec1, const std::array<scalar, 2>& vec2) {
    scalar diffX = vec1[0] - vec2[0];
    scalar diffY = vec1[1] - vec2[1];
    return std::sqrt(diffX * diffX + diffY * diffY);
}

// Function representing the exact solution
scalar uExact(const std::array<scalar, 2>& x) {
    std::array<scalar, 2> a = {-16.0 / 15.0, 0.0};
    std::array<scalar, 2> b = {-1.0 / 15.0, 0.0};

    return (std::log(norm(x, a)) - std::log(norm(x, b))) / std::log(2) - 1;
}

scalar L2_norm(const Mesh &mesh, scalar *sol) {
    scalar errorSum = 0.0;
    for (int i = 0; i < mesh.n_nodes; ++i) {
        // Assuming each entry in 'approx' corresponds to the x-coordinate of a point on the x-axis
        std::array<scalar, 2> vertex{mesh.vertices.at(i).x, mesh.vertices.at(i).y};// Construct the point (x, 0)
        scalar exactValue = uExact(vertex); // Evaluate the exact solution
        scalar error = exactValue - sol[i]; // Compute the difference
        errorSum += error * error; // Add the square of the difference to the sum
    }
    return sqrt(errorSum);
}


int main() {

    std::chrono::time_point< std::chrono::high_resolution_clock > t_start, t_end;
    
    std::ofstream output("time_series.txt");
    assert(output.is_open(),
            "Failed to open the mesh file");
    
    std::string file_name = "n_1h_circle";
    std::string folder_name = "1h_linear_circle_series/";

    for(int i = 4; i <= MAX_N; i += 2) {
        //t_start = std::chrono::high_resolution_clock::now();
        auto mesh = Mesh::parse_mesh("../meshes/" + folder_name + std::to_string(i) + file_name + ".msh");
        auto res = solve_fem(mesh);
        //t_end = std::chrono::high_resolution_clock::now();
        //output << static_cast<std::chrono::duration<scalar> >(t_end - t_start).count() << ", ";
        output << L2_norm(mesh, res.data) << ", ";
    }


    return 0;
}
