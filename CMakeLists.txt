cmake_minimum_required(VERSION 3.5)
project(ico)

# ------------------
#   Compiler Flags
# ------------------

set (CMAKE_CXX_FLAGS " -O2 -march=native -std=c++17 -Wno-unused-parameter -Wno-sign-compare ")

#add_subdirectory(xdrfile-1.1.4)


add_executable (ico main.cpp 
xdrfile-1.1.4/include/xdrfile.h 
xdrfile-1.1.4/include/xdrfile_trr.h 
xdrfile-1.1.4/include/xdrfile_xtc.h 
xdrfile-1.1.4/src/xdrfile.c 
xdrfile-1.1.4/src/xdrfile_c_test.c 
xdrfile-1.1.4/src/xdrfile_trr.c 
xdrfile-1.1.4/src/xdrfile_xtc.c
particle.h 
icosahedron.h 
oblatespheroid.h 
sphere.h 
tennisball.h 
spherepatch.h 
pentamer.h 
dodecahedron.h 
surface.h 
chain.h 
data.h 
slab.h 
atom.h 
cow.h 
rng.h 
atom.h 
chain.h 
cow.h 
data.h 
dodecahedron.h 
icosahedron.h 
oblatespheroid.h 
particle.h 
pentamer.h 
rng.h 
slab.h 
sphere.h 
spherepatch.h 
surface.h 
tennisball.h 
xtcanalysis.h 
welford.h
force_field.h
)

target_include_directories(ico PUBLIC "xdrfile-1.1.4/include/")

