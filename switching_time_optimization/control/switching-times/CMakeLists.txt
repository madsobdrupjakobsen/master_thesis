# cmake_minimum_required(VERSION 3.10) set(CMAKE_C_COMPILER
# /usr/local/bin/gcc-9) set(CMAKE_CXX_COMPILER /usr/local/bin/g++-9)
# project(SwitchingTimes) set(CMAKE_CXX_STANDARD 17) set(IPOPT_LIBRARY_DIRS
# /usr/local/lib) set(IPOPT_INCLUDE_DIRS /usr/local/include/coin-or)
# include_directories("/usr/local/include")
# include_directories(${IPOPT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
# link_directories(${IPOPT_LIBRARY_DIRS}) find_package(OpenMP REQUIRED)
# add_executable(SwitchingTimes main.cpp src/switching-times.hpp src/switching-
# times.cpp src/cppad-eigen.hpp src/cppad-eigen-odeint.hpp src/switching-times-
# example.cpp)

# target_link_libraries(SwitchingTimes PRIVATE OpenMP::OpenMP_CXX)
# target_link_libraries(SwitchingTimes PRIVATE ipopt)

# BELOW is pybind11 CMakeLists.txt - above is ordinary CMakeLists.txt

cmake_minimum_required(VERSION 3.10)
set(CMAKE_C_COMPILER /usr/local/bin/gcc-9)
set(CMAKE_®CXX_COMPILER /usr/local/bin/g++-9)
project(switching_times)
set(CMAKE_CXX_STANDARD 17)
set(IPOPT_LIBRARY_DIRS /usr/local/lib)
set(IPOPT_INCLUDE_DIRS /usr/local/include/coin-or)
include_directories("/usr/local/include")
include_directories(${IPOPT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
link_directories(${IPOPT_LIBRARY_DIRS})
include_directories("./pybind11/include")
add_subdirectory(pybind11)

# ------- OG modules
# First order model
pybind11_add_module(
  switching_times_1st
  main_1st.cpp
  src/switching-times.hpp
  src/switching-times.cpp
  src/cppad-eigen.hpp
  src/cppad-eigen-odeint.hpp
  src/switching-times-ice-tank-simple.cpp)
target_link_libraries(switching_times_1st PRIVATE ipopt)

# Second order model
pybind11_add_module(
  switching_times_2nd
  main_2nd_order.cpp
  src/switching-times.hpp
  src/switching-times.cpp
  src/cppad-eigen.hpp
  src/cppad-eigen-odeint.hpp
  src/switching-times-ice-tank-2ndorder.cpp)
target_link_libraries(switching_times_2nd PRIVATE ipopt)

# Third order model
pybind11_add_module(
  switching_times_3rd
  main_3rd_order.cpp
  src/switching-times.hpp
  src/switching-times.cpp
  src/cppad-eigen.hpp
  src/cppad-eigen-odeint.hpp
  src/switching-times-ice-tank-3rdorder.cpp)
target_link_libraries(switching_times_3rd PRIVATE ipopt)

# ------- RK modules

# First order model
pybind11_add_module(
  switching_times_1st_rk
  main_1st_rk.cpp
  src/switching-times-rk.hpp
  src/switching-times-rk.cpp
  src/cppad-eigen.hpp
  src/cppad-eigen-odeint.hpp
  src/switching-times-ice-tank-simple-rk.cpp)
target_link_libraries(switching_times_1st_rk PRIVATE ipopt)

# Second order model
pybind11_add_module(
  switching_times_2nd_rk
  main_2nd_rk.cpp
  src/switching-times-rk.hpp
  src/switching-times-rk.cpp
  src/cppad-eigen.hpp
  src/cppad-eigen-odeint.hpp
  src/switching-times-ice-tank-2ndorder_rk.cpp)
target_link_libraries(switching_times_2nd_rk PRIVATE ipopt)

# Third order model
pybind11_add_module(
  switching_times_3rd_rk
  main_3rd_rk.cpp
  src/switching-times-rk.hpp
  src/switching-times-rk.cpp
  src/cppad-eigen.hpp
  src/cppad-eigen-odeint.hpp
  src/switching-times-ice-tank-3rdorder_rk.cpp)
target_link_libraries(switching_times_3rd_rk PRIVATE ipopt)

