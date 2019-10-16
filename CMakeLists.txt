cmake_minimum_required(VERSION 3.14)
project(cse570project2)

set(CMAKE_CXX_STANDARD 17)

include_directories(./Core)
#include_directories(/Library/Developer/CommandLineTools/usr/include/c++/v1)


add_library(halfedge
        ./Core/Edge.cpp
        ./Core/Face.cpp
        ./Core/FException.cpp
        ./Core/HalfEdge.cpp
        ./Core/OBJFileReader.cpp
        ./Core/Parser.cpp
        ./Core/Point.cpp
        ./Core/Solid.cpp
        ./Core/SolidDelegate.cpp
        ./Core/StringTokenizer.cpp
        ./Core/TopologyException.cpp
        ./Core/Trait.cpp
        ./Core/Vertex.cpp)

add_executable(cse570project2 Subdivision/Subdivision.cpp)

target_link_libraries(cse570project2 halfedge)