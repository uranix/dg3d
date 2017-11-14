cmake_minimum_required(VERSION 3.0)

project(dg3d)

add_subdirectory(libMeshb)

add_executable(gmf2vtk gmf2vtk.cpp)
include_directories(libMeshb/sources/)
target_link_libraries(gmf2vtk Meshb.7)

add_executable(main main.cpp mesh.cpp)
include_directories(libMeshb/sources/)
target_link_libraries(main Meshb.7)
