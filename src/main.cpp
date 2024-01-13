#include "FEM.h"

#include "benchmark.h"

#define PI 3.14159

scalar source_fn(const Vector<2> &x) {
    return (8.0 * PI * PI + 1) 
        * std::cos(2 * PI * x(0)) 
        * std::cos(2 * PI * x(1));
}


int main() {


    auto mesh = Mesh::load("../meshes/first_mesh_n100.msh");
    auto res = solve_fem(mesh, &source_fn);


    benchmark::global_timer::write_to_file("benchmark.json");

    return 0;
}

