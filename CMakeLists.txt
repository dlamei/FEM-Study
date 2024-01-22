cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(FEM-Study)

add_subdirectory(libs/eigen)

set(SOURCES
	src/utils.h
    src/benchmark.h
		src/geometry.h
		src/geometry.cpp

    src/FEM.h src/FEM.cpp
	src/FEM_conj_grad.h
	src/FEM_conj_grad.cpp
	
	src/MatrixImpl/conj_grad.h 
	src/MatrixImpl/conj_grad.cpp
	src/MatrixImpl/conj_grad.ocl
	
    src/MatrixImpl/stack_matrix.h
    src/MatrixImpl/stack_matrix.cpp
    src/MatrixImpl/naive_matrix.h
    src/MatrixImpl/naive_matrix.cpp
    src/MatrixImpl/band_matrix.h
    src/MatrixImpl/band_matrix.cpp
    src/MatrixImpl/sparse_matrix.h
    src/MatrixImpl/sparse_matrix.cpp
	)

set(SOURCES_NAIVE
	src/utils.h
    src/benchmark.h
		src/geometry.h
		src/geometry.cpp
	src/FEM_naive.h
	src/FEM_naive.cpp

    src/MatrixImpl/stack_matrix.h
    src/MatrixImpl/stack_matrix.cpp
    src/MatrixImpl/naive_matrix.h
    src/MatrixImpl/naive_matrix.cpp
	)

find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})

add_executable(FEM ${SOURCES} src/main.cpp)
add_executable(naive_FEM ${SOURCES_NAIVE} src/main_naive.cpp)
add_executable(time_series ${SOURCES_NAIVE} src/time_series.cpp)

# Add optimization flags to the target
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(time_series PRIVATE -O3)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(time_series PRIVATE /O2)
endif()
# Add optimization flags to the target
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(naive_FEM PRIVATE -O3)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(naive_FEM PRIVATE /O2)
endif()


target_link_libraries(FEM eigen)

target_link_libraries(time_series eigen)

target_link_libraries(FEM ${OpenCL_LIBRARIES})



# testing
add_executable(FEM_TEST EXCLUDE_FROM_ALL ${SOURCES} src/test.cpp)
target_compile_definitions(FEM_TEST PRIVATE COMPILE_TESTS)
find_package(Threads)
target_link_libraries(FEM_TEST ${CMAKE_THREAD_LIBS_INIT} eigen)

target_link_libraries(FEM_TEST ${OpenCL_LIBRARIES})

# profiling
add_executable(FEM_PROFILE ${SOURCES} src/main.cpp)
target_link_libraries(FEM_PROFILE eigen)
target_compile_options(FEM_PROFILE PRIVATE -pg)
target_link_options(FEM_PROFILE PRIVATE -pg)

# benchmark
add_executable(FEM_BENCH ${SOURCES} src/main.cpp)
target_link_libraries(FEM_BENCH eigen)
target_compile_definitions(FEM_BENCH PRIVATE PROFILING=1)


# custom targets 
add_custom_target(r DEPENDS run)
add_custom_target(run
	COMMAND echo "running FEM..." && echo "" && ./FEM
	DEPENDS FEM
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	)

add_custom_target(t DEPENDS test)
add_custom_target(test
	COMMAND echo "running tests..." && echo "" && ./FEM_TEST
	DEPENDS FEM_TEST
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	)

add_custom_target(p DEPENDS profile)
add_custom_target(profile
	COMMAND echo "running profiling..." && echo "" && ./FEM_PROFILE
	COMMAND ../scripts/run_profiling.sh
	DEPENDS FEM_PROFILE
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	)

add_custom_target(bench DEPENDS benchmark)
add_custom_target(benchmark
    COMMAND echo "running benchmark..." && echo "" && ./FEM_BENCH
    DEPENDS FEM_BENCH
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
	)
