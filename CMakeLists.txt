cmake_minimum_required(VERSION 3.10.2)
project(EMplanner)
set(CMAKE_BUILD_TYPE "Debug")

list(APPEND CMAKE_PREFIX_PATH "/home/linkx/local_install_packages/foxglove-websocket")
list(APPEND CMAKE_PREFIX_PATH "/home/linkx/local_install_packages/flatbuffers")
list(APPEND CMAKE_PREFIX_PATH "/home/linkx/local_install_packages/nlohmann_json")
list(APPEND CMAKE_PREFIX_PATH "/home/linkx/local_install_packages/websocketpp")

find_package(Eigen3 REQUIRED)
find_package(osqp REQUIRED)
find_package(OsqpEigen REQUIRED)

find_package(foxglove-schemas-protobuf CONFIG REQUIRED)
find_package(foxglove-websocket REQUIRED)
find_package(nlohmann_json REQUIRED)
find_package(websocketpp REQUIRED)
find_package(Boost REQUIRED COMPONENTS system)
find_package(OpenSSL REQUIRED)
find_package(ZLIB REQUIRED)



include_directories(/usr/local/include/coin-or ./include)
link_directories(/usr/local/lib ./lib)

add_executable(main src/main.cpp 
                    lib/smoothSolver.cpp 
                    lib/frenetUtils.cpp 
                    lib/osqpSolver.cpp 
                    lib/server_config.cpp
                    lib/visualization.cpp)

target_link_libraries(main PUBLIC 
                      ipopt 
                      coinhsl 
                      coinmumps 
                      lapack 
                      blas 
                      Eigen3::Eigen 
                      osqp::osqp 
                      OsqpEigen::OsqpEigen)

target_link_libraries(main PUBLIC
    foxglove-schemas-protobuf::foxglove-schemas-protobuf
    foxglove-websocket::foxglove-websocket
)

target_link_libraries(main PUBLIC 
                      nlohmann_json::nlohmann_json 
                      ZLIB::ZLIB
                      websocketpp::websocketpp 
                      Boost::system 
                      pthread
                      OpenSSL::SSL 
                      OpenSSL::Crypto)

