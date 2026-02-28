cmake_minimum_required(VERSION 3.10)
project(tksvm-client CXX)

add_compile_options(
  "-Wall" "-Wextra"
  "-Wfatal-errors"
  "-Wno-unknown-pragmas"
  "-Wno-c++17-extensions"
  "-Wno-unused-parameter"
  "-fdiagnostics-color=always"
  )

find_package(ALPSCore 2.2.0 REQUIRED)
find_package(Threads)

add_subdirectory(tk-svm)

include_directories(include ${TKSVM_INCLUDE_DIRS})

add_compile_definitions(
  ${TKSVM_DEFINITIONS}
  TKSVM_SIMINCL=client/sim.hpp
  )

#Add all relevant cpp files of the client code
#set(CLIENT_SRC
#    )

add_executable(sample ${TKSVM_SAMPLE_SRC} ${CLIENT_SRC})
add_executable(learn ${TKSVM_LEARN_SRC} ${CLIENT_SRC})
#add_executable(test ${TKSVM_TEST_SRC} ${CLIENT_SRC})
add_executable(coeffs ${TKSVM_COEFFS_SRC} ${CLIENT_SRC})
add_executable(segregate-phases ${TKSVM_SEGREGATE_PHASES_SRC} ${CLIENT_SRC})

target_link_libraries(sample ${ALPSCore_LIBRARIES} ${TKSVM_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(learn ${ALPSCore_LIBRARIES} ${TKSVM_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#target_link_libraries(test ${ALPSCore_LIBRARIES} ${TKSVM_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(coeffs ${ALPSCore_LIBRARIES} ${TKSVM_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(segregate-phases ${ALPSCore_LIBRARIES} ${TKSVM_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS
  sample
  learn
  #test
  coeffs
  segregate-phases
  DESTINATION bin)

