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
cmake .. -DPYTHON_EXECUTABLE=$(which python) -DCMAKE_BUILD_TYPE=Release