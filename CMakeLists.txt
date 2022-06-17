cmake_minimum_required(VERSION 3.10.2)
project(EMplanner)
set(CMAKE_BUILD_TYPE "Debug")


find_package(Eigen3 REQUIRED)
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Python3 COMPONENTS NumPy)


include_directories(/usr/local/include/coin-or ./include)
link_directories(/usr/local/lib ./lib)

add_executable(main src/main.cpp lib/smoothSolver.cpp lib/frenetUtils.cpp lib/osqpSolver.cpp)

target_link_libraries(main ipopt coinhsl coinmumps lapack blas Eigen3::Eigen osqp::osqp OsqpEigen::OsqpEigen   Python3::Python
         Python3::Module Python3::NumPy
 )
