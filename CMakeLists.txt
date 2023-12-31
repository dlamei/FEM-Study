cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(FEM-Study)

add_subdirectory(libs/eigen)

set(SOURCES
	src/utils.h
	src/Geometries/geometry.h 
	src/Geometries/geometry.cpp
	src/matrix.cpp src/matrix.h
    src/FEM.h src/FEM.cpp
	)


add_executable(FEM ${SOURCES} src/main.cpp)
target_link_libraries(FEM eigen)

# testing
add_executable(FEM_TEST EXCLUDE_FROM_ALL ${SOURCES} src/test.cpp)
target_compile_definitions(FEM_TEST PRIVATE COMPILE_TESTS)
find_package(Threads)
target_link_libraries(FEM_TEST ${CMAKE_THREAD_LIBS_INIT} eigen)



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
