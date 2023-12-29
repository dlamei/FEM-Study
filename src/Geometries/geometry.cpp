#include "geometry.h"

using namespace Geometry;

Mesh Mesh::parse_mesh(std::string file_name) {
    Mesh mesh{};

    std::fstream input;
    input.open(file_name);
    assert(input.is_open(), 
            "Failed to open the mesh file");
    
    // read in number of vertices, triangeles and edge vertices
    input >> mesh.nof_vertices;
    input >> mesh.nof_triangles;
    input >> mesh.nof_boundry_edges;
    assert((mesh.nof_vertices >= 0 && mesh.nof_triangles >= 0 && mesh.nof_boundry_edges >= 0), 
            "Error nof_vertices, nof_triangles or nof_boundry_edges is an invalid input");
    
    // read in vertice coordinates
    for(int i = 0; i < mesh.nof_vertices; ++i) {
        scalar x, y;
        u64 boundry_label;
        input >> x;
        input >> y;
        input >> boundry_label;
        mesh.vertices.push_back({x,y});
    }
    // read in the triangles
    for(int i = 0; i < mesh.nof_triangles; ++i) {
        u64 l, m, n;
        u64 boundry_label;
        input >> l;
        input >> m;
        input >> n;
        input >> boundry_label;
        mesh.triangles.push_back({l - 1, m - 1, n - 1});
    }
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
    std::cout << "Vertices: " << mesh.nof_vertices << "   ";
    std::cout << "Triangles: " << mesh.nof_triangles << "   ";
    std::cout << "Boundry: " << mesh.nof_boundry_edges << '\n';
    return mesh;
}

matrix_t Mesh::get_tria_coords(int i) {
    matrix_t tria_coords;
    tria_coords.push_back(vertices.at(triangles.at(i).at(0)));
    tria_coords.push_back(vertices.at(triangles.at(i).at(1)));
    tria_coords.push_back(vertices.at(triangles.at(i).at(2)));
    return tria_coords;
}

void Mesh::save_mesh_3D(std::string filename, const std::vector<scalar> &z) const {
    std::ofstream output(filename);
    assert(output.is_open(),
            "Failed to open the mesh file");
    output << "# vtk DataFile Version 3.0\n";
    output << "2D PDE Solution on Triangulated Mesh\n";
    output << "ASCII\n";
    output << "DATASET UNSTRUCTURED_GRID\n";
    output << "POINTS " << nof_vertices << " float\n";
            
    for (int i = 0; i < nof_vertices; ++i) {
        output << vertices.at(i).at(0) << ' ' << vertices.at(i).at(1) << ' ' << "0" << '\n';
    }
    output << '\n';
    output << "CELLS " << nof_triangles << " " << nof_triangles * 4 << "\n";
    for (int i = 0; i < nof_triangles; ++i) {
        output << "3 " << triangles.at(i).at(0) << " " << triangles.at(i).at(1) << " " << triangles.at(i).at(2) << '\n'; 
    }
    output << '\n';
    output << "CELL_TYPES " << nof_triangles << '\n';
    for (int i = 0; i < nof_triangles; ++i) {
        output << "5\n";
    }
    output << '\n';
    output << "POINT_DATA " << nof_vertices << '\n';
    output << "SCALARS solution_field float 1\n";
    output << "LOOKUP_TABLE default\n";
    for (int i = 0; i < nof_vertices; ++i) {
        output << z.at(i) << '\n';
    }
    // Print message
    std::cout << "Finished saving Solution to " + filename << '\n';
    std::cout << "Vertices: " << nof_vertices << "   ";
    std::cout << "Triangles: " << nof_triangles << "   ";
    std::cout << "Boundry: " << nof_boundry_edges << '\n';
}
