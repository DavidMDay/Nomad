
cmake_minimum_required(VERSION 3.25)
set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall")
set(CMAKE_CXX_FLAGS_RELEASE_INIT "-Wall")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
project(matrix_computations)


# Add library
add_library(matrixlib util.cpp QRdecomposition.cpp hqr2.cpp applyQ.cpp applyQT.cpp triuSol.cpp)

set_property(TARGET matrixlib PROPERTY CXX_STANDARD 17)

# Add executable
add_executable(tester main.cpp stridedInnerProductTest.cpp qrTest.cpp matrixCollection.cpp testDecomposition.cpp getSquareUpperTri.cpp triSolTest.cpp)

# Link executable to library
target_link_libraries(tester matrixlib)

target_compile_options(tester PUBLIC -std=c++14 -Wall)
set_property(TARGET tester PROPERTY CXX_STANDARD 17)

