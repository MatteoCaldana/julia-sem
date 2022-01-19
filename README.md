# deal-rl



## DealLES
export mkPrefix=/u/sw
source $mkPrefix/etc/profile
module load gcc-glibc/11 dealii

cd dealii/
mkdir build
cd build/
cmake -DCMAKE_BUILD_TYPE="Release" .. 