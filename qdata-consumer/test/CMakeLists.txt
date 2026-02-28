cmake_minimum_required(VERSION 3.10)
project(lattice_test CXX)

SET( CMAKE_lattice_test_COMPILER "g++" )

add_compile_options(
  "-Wall" "-Wextra"
  "-Wfatal-errors"
  "-Wno-unknown-pragmas"
  "-Wno-c++17-extensions"
  "-fdiagnostics-color=always"
  )

find_package(ALPSCore 2.2.0 REQUIRED)
find_package(Threads)

add_subdirectory(../tk-svm tksvm)

include_directories(../include ${TKSVM_INCLUDE_DIRS})
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

add_executable(test_lattice lattice.cpp)

target_link_libraries(test_lattice ${ALPSCore_LIBRARIES} ${TKSVM_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS
    test_lattice
  DESTINATION bin)

