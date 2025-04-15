# ./Helper.sh
# export PATH="$PATH:/Applications/CMake.app/Contents/bin"
# cmake -DCMAKE_BUILD_TYPE=Debug  ../Nomad/matrix_computations
# mv compile_commands.json ../Nomad/matrix_computations
# cmake --build .
# ./tester
# ctest
# mv libmatrixlib.a ../matrix_computations
 rm libmatrixlib.a compile_commands.json
 rm tester first CMakeCache.txt Makefile cmake_install.cmake
 rm CTestTestfile.cmake
 rm -rf CMakeFiles Testing 
