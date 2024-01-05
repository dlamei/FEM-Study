#include "geometry_eigen.h"

using namespace Geometry;

Mesh Mesh::parse_mesh(std::string file_name) {
    Mesh mesh{};

    std::fstream input;
    input.open(file_name);
    assert(input.is_open(), 
            "Failed to open the mesh file");

    // read in number of nodes, triangeles and edge nodes
    input >> mesh.nof_nodes;
    input >> mesh.nof_triangles;
    input >> mesh.nof_boundry_edges;
    assert((mesh.nof_nodes >= 0 && mesh.nof_triangles >= 0 && mesh.nof_boundry_edges >= 0), 
            "Error nof_nodes, nof_triangles or nof_boundry_edges is an invalid input");

    // read in node coordinates
    auto nodes = Eigen::Matrix<scalar, Eigen::Dynamic, 2>(mesh.nof_nodes, 2);
    for(int i = 0; i < mesh.nof_nodes; ++i) {
        scalar x, y;
        u64 boundry_label;
        input >> x;
        input >> y;
        input >> boundry_label;
        nodes(i, 0) = x;
        nodes(i, 1) = y;
    }
    mesh.nodes = nodes;

    // read in the triangles
    auto triangles = Eigen::Matrix<int, Eigen::Dynamic, 3>(mesh.nof_triangles, 3);
    for(int i = 0; i < mesh.nof_triangles; ++i) {
        u64 l, m, n;
        u64 boundry_label;
        input >> l;
        input >> m;
        input >> n;
        input >> boundry_label;
        triangles(i, 0) = l - 1;
        triangles(i, 1) = m - 1;
        triangles(i, 2) = n - 1;
    }
    mesh.triangles = triangles;

    // read in the boundries
    mesh.boundries.push_back( {} );
    for(int i = 0; i < mesh.nof_boundry_edges; ++i) {
        u64 m, n;
        u64 boundry_label;
        input >> m;
        input >> n;
        input >> boundry_label;
        mesh.boundries.at(0).push_back( {m - 1, n - 1} );
        while( boundry_label >= mesh.boundries.size() ) {
            mesh.boundries.push_back( {} );
        }
        mesh.boundries.at(boundry_label).push_back( {m - 1, n - 1} );
    }

    // Print message
    std::cout << "Finished parsing Mesh from " + file_name + "\n";
    std::cout << "nodes: " << mesh.nof_nodes << "   ";
    std::cout << "Triangles: " << mesh.nof_triangles << "   ";
    std::cout << "Boundry: " << mesh.nof_boundry_edges << '\n';
    return mesh;
}


Triangle Mesh::get_tria_coords(int i) const {
    Triangle tria_coords;
    tria_coords << nodes(triangles(i, 0), 0),
                   nodes(triangles(i, 1), 0),
                   nodes(triangles(i, 2), 0),
                   nodes(triangles(i, 0), 1),
                   nodes(triangles(i, 1), 1),
                   nodes(triangles(i, 2), 1);
    return tria_coords;
}

void Mesh::save_mesh_3D(std::string filename, const Vector<Eigen::Dynamic> &z) const {
    std::ofstream output(filename);
    assert(output.is_open(),
            "Failed to open the mesh file");
    output << "# vtk DataFile Version 3.0\n";
    output << "2D PDE Solution on Triangulated Mesh\n";
    output << "ASCII\n";
    output << "DATASET UNSTRUCTURED_GRID\n";
    output << "POINTS " << nof_nodes << " float\n";

    for (int i = 0; i < nof_nodes; ++i) {
        output << nodes(i, 0) << ' ' << nodes(i, 1) << ' ' << "0" << '\n';
    }
    output << '\n';
    output << "CELLS " << nof_triangles << " " << nof_triangles * 4 << "\n";
    for (int i = 0; i < nof_triangles; ++i) {
        output << "3 " << triangles(i, 0) << " " << triangles(i, 1) << " " << triangles(i, 2) << '\n'; 
    }
    output << '\n';
    output << "CELL_TYPES " << nof_triangles << '\n';
    for (int i = 0; i < nof_triangles; ++i) {
        output << "5\n";
    }
    output << '\n';
    output << "POINT_DATA " << nof_nodes << '\n';
    output << "SCALARS solution_field float 1\n";
    output << "LOOKUP_TABLE default\n";
    for (int i = 0; i < nof_nodes; ++i) {
        output << z(i)  << '\n';
    }
    // Print message
    std::cout << "Finished saving Solution to " + filename << '\n';
    std::cout << "nodes: " << nof_nodes << "   ";
    std::cout << "Triangles: " << nof_triangles << "   ";
    std::cout << "Boundry: " << nof_boundry_edges << '\n';
}
