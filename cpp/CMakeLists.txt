cmake_minimum_required(VERSION 3.9)
project(heat)

set(CMAKE_CXX_STANDARD 11)

include_directories(/opt/ipopt/include/coin)    # Ipopt header files are here
link_directories(/opt/ipopt/lib                 # Ipopt libraries
        /home/bt304019/Downloads/OpenBLAS-0.2.19)   # My version of Ipopt requires OpenBLAS which is here

add_executable(heat
        heat_nlp.cpp
        heat_nlp.hpp
        main.cpp
        matrixop.cpp
        matrixop.hpp)

target_link_libraries(heat ipopt boost_system boost_filesystem)
