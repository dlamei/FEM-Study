# FEM-Study — Methods

Summary of implemented numerical methods and where to find them.

- FEM formulation
  - Galerkin finite element method (2D Poisson-style problems)
  - Linear (P1 / barycentric) triangular elements — element stiffness & local load-vector assembly
  - Assembly routines: src/FEM_naive.cpp, src/FEM.cpp, src/FEM_conj_grad.cpp

- Boundary conditions & I/O
  - Dirichlet boundary application (apply_boundry_conditions)
  - Mesh I/O: reads Gmsh .msh, writes VTK for visualization
  - See: src/FEM_naive.cpp (apply_boundry_conditions), geometry/mesh utilities

- Matrix backends
  - Naive dense/custom matrix (naive_matrix)
  - Band, stack, and sparse custom matrix implementations
  - Locations: src/MatrixImpl/{naive_matrix,band_matrix,stack_matrix,sparse_matrix}.*

- Linear solvers
  - Direct:
    - Custom LDLT — used by naive backend
    - Eigen sparse direct solvers (SimplicialLDLT, SparseLU) — src/FEM.cpp
  - Iterative:
    - Conjugate Gradient (CPU) — src/MatrixImpl/conj_grad.*
    - OpenCL-accelerated Conjugate Gradient (fem_ocl::run_conj_grad) — src/MatrixImpl/conj_grad.ocl, used by src/FEM_conj_grad.cpp

- Tooling hints (methods-level)
  - Solver selection is compile-time in src/main.cpp (FEM_EIGEN_IMPL, FEM_CONJ_GRAD_IMPL, FEM_NAIVE_IMPL)
  - Simple benchmarking / time-series harness available: src/time_series.cpp

This project is a compact playground for comparing assembly strategies, matrix storage formats, and solver approaches (direct vs iterative, CPU vs OpenCL).
