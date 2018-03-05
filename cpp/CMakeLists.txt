cmake_minimum_required(VERSION 3.8)
project(heat)

set(CMAKE_CXX_STANDARD 11)

include_directories(/opt/ipopt/include/coin # Ipopt header files are here
                    /opt/OpenBLAS/include)
link_directories(/opt/ipopt/lib                 # Ipopt libraries
                 /opt/OpenBLAS/lib)   # My version of Ipopt requires OpenBLAS which is here

add_executable(heat
        heat_nlp.cpp
        heat_nlp.hpp
        main.cpp
        matrixop.cpp
        matrixop.hpp)

target_link_libraries(heat ipopt lapack blas boost_system boost_filesystem)
