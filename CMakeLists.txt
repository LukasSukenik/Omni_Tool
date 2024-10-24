cmake_minimum_required(VERSION 3.5)
project(Omni_Tool)

# ------------------
#   Compiler Flags
# ------------------

set (CMAKE_CXX_FLAGS " -O2 -march=native -std=c++17 -Wno-unused-parameter -Wno-sign-compare ")

#add_subdirectory(xdrfile-1.1.4)


add_executable (Omni_Tool main.cpp 
xdrfile-1.1.4/include/xdrfile.h 
xdrfile-1.1.4/include/xdrfile_trr.h 
xdrfile-1.1.4/include/xdrfile_xtc.h 
xdrfile-1.1.4/src/xdrfile.c 
xdrfile-1.1.4/src/xdrfile_c_test.c 
xdrfile-1.1.4/src/xdrfile_trr.c 
xdrfile-1.1.4/src/xdrfile_xtc.c
particles/chain.h
particles/cow.h
particles/dodecahedron.h
particles/ellipsoid.h
particles/globular_sphere.h
particles/icosahedron.h
particles/oblatespheroid.h
particles/pentamer.h
particles/sphere.h
particles/spherepatch.h
particles/surface.h
particles/tennisball.h
IO/io_input.h
IO/io_lammps.h
IO/io_pdb.h
angle.h
atom.h
bond.h
data.h
force_field.h
particle.h
rng.h
sim_box.h
slab.h
types.h
xtcanalysis.h
welford.h
)

target_include_directories(Omni_Tool PUBLIC "xdrfile-1.1.4/include/" "particles/" "IO/")

