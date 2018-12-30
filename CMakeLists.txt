cmake_minimum_required(VERSION 3.13)
project(pcb_sat)

set(CMAKE_CXX_STANDARD 17)


FIND_PACKAGE(cryptominisat5)
FIND_PACKAGE(Boost COMPONENTS system iostreams filesystem)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)



INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
add_executable(pcb_sat main.cpp Expression.cpp Expression.h Solg_solver.cpp Solg_solver.h Solg.cpp Solg.h Circuit.cpp Circuit.h Dynamic.cpp Dynamic.h)
target_compile_options(pcb_sat PUBLIC -O3 -march=native)
target_link_libraries(pcb_sat ${CRYPTOMINISAT5_LIBRARIES} ${Boost_LIBRARIES} Eigen3::Eigen xxhash)