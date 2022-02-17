# deal-rl



## DealLES
export mkPrefix=/u/sw
source $mkPrefix/etc/profile
module load gcc-glibc/11 dealii pybind11

cd dealii/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE="Release" .. 
cmake --build . --config Release

## myfex
cmake .. -DPYTHON_EXECUTABLE=/u/sw/toolchains/gcc-glibc/11.2.0/base/bin/python -DCMAKE_BUILD_TYPE=Release