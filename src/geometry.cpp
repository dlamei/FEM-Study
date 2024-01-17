#include "geometry.h"

Mesh Mesh::parse_mesh(std::string file_name) {
    Mesh mesh{};

    std::fstream input;
    input.open(file_name);
    assert(input.is_open(), 
            "Failed to open the mesh file");

    // read in number of vertices, triangeles and edge vertices
    input >> mesh.n_nodes;
    input >> mesh.n_triangles;
    input >> mesh.n_boundry_nodes;
    assert((mesh.n_nodes >= 0 && mesh.n_triangles >= 0 && mesh.n_boundry_nodes >= 0), 
            "Error n_nodes, n_triangles or n_boundry_nodes is an invalid input");

    // read in vertice coordinates
    for(int i = 0; i < mesh.n_nodes; ++i) {
        scalar x, y;
        index_t boundry_label;
        input >> x;
        input >> y;
        input >> boundry_label;
        mesh.vertices.push_back({x,y});
    }
    // read in the triangles
    for(int i = 0; i < mesh.n_triangles; ++i) {
        index_t l, m, n;
        index_t boundry_label;
        input >> l;
        input >> m;
        input >> n;
        input >> boundry_label;
        mesh.triangles.push_back({l - 1, m - 1, n - 1});
    }
    // read in the boundries
    for(int i = 0; i < mesh.n_boundry_nodes; ++i) {
        index_t m, n;
        index_t boundry_label;
        input >> m;
        input >> n;
        input >> boundry_label;
        switch( boundry_label ) {
            case 1: 
                mesh.outer_boundries.push_back( m - 1 );
                break;
            case 2:
                mesh.inner_boundries.push_back( m - 1 );
                break;
            default:
                assert(false, "Error, too many different boundries, change switch statement in parse_mesh");

        }
    }

    // Print message
    std::cout << "Finished parsing Geometry from " + file_name + "\n";
    std::cout << "Vertices: " << mesh.n_nodes << "   ";
    std::cout << "Triangles: " << mesh.n_triangles << "   ";
    std::cout << "Boundry: " << mesh.n_boundry_nodes << '\n';
    return mesh;
}

void Mesh::save_mesh_3D(const std::string& filename, const scalar *result, usize count) const {
    assert(count == n_nodes, "sanity check");

    std::ofstream output(filename);
    assert(output.is_open(),
            "Failed to open the mesh file");
    output << "# vtk DataFile Version 3.0\n";
    output << "2D PDE Solution on Triangulated Geometry\n";
    output << "ASCII\n";
    output << "DATASET UNSTRUCTURED_GRID\n";
    output << "POINTS " << n_nodes << " float\n";

    for (int i = 0; i < n_nodes; ++i) {
        output << vertices.at(i).x << ' ' << vertices.at(i).y << ' ' << "0" << '\n';
    }
    output << '\n';
    output << "CELLS " << n_triangles << " " << n_triangles * 4 << "\n";
    for (int i = 0; i < n_triangles; ++i) {
        output << "3 " << triangles.at(i).a << " " << triangles.at(i).b << " " << triangles.at(i).c << '\n'; 
    }
    output << '\n';
    output << "CELL_TYPES " << n_triangles << '\n';
    for (int i = 0; i < n_triangles; ++i) {
        output << "5\n";
    }
    output << '\n';
    output << "POINT_DATA " << n_nodes << '\n';
    output << "SCALARS solution_field float 1\n";
    output << "LOOKUP_TABLE default\n";
    for (int i = 0; i < n_nodes; ++i) {
        output << result[i] << '\n';
    }
    // Print message
    std::cout << "Finished saving Solution to " + filename << '\n';
    std::cout << "Vertices: " << n_nodes << "   ";
    std::cout << "Triangles: " << n_triangles << "   ";
    std::cout << "Boundry: " << n_boundry_nodes << '\n';
}
