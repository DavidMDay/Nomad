cmake_minimum_required(VERSION 3.24.1)
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Wall")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
project(harminc_balance_project)

# Add executable
add_executable(tester main.cpp Tests.cpp FourierAnalysis.cpp IntegratorTest.cpp TestOperators.cpp)
# Add library
add_library(hblib FourierCoefficients.cpp SpaceTimeTransient.cpp FixedStepIntegrator.cpp HB.cpp Operators.cpp matrix2x2.cpp smp.cpp)

set_property(TARGET hblib PROPERTY CXX_STANDARD 17)

add_library(libmatrixlib STATIC IMPORTED)
include_directories(libmatrixlib INTERFACE ../matrix_computations)
set_target_properties(libmatrixlib PROPERTIES IMPORTED_LOCATION ../matrix_computations/libmatrixlib.a)

target_link_libraries(tester hblib libmatrixlib)

target_compile_options(tester PUBLIC -std=c++14 -Wall)
set_property(TARGET tester PROPERTY CXX_STANDARD 17)

