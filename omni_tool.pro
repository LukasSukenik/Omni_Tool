TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += xdrfile-1.1.4/include/ particles/ IO/ Systems/

QMAKE_CXXFLAGS += -O2 -march=native -std=c++17 -Wno-unused-parameter -Wno-sign-compare

SOURCES += main.cpp \
    xdrfile-1.1.4/src/xdrfile.c \
    xdrfile-1.1.4/src/xdrfile_c_test.c \
    xdrfile-1.1.4/src/xdrfile_trr.c \
    xdrfile-1.1.4/src/xdrfile_xtc.c

HEADERS += \
    angle.h \
    atom.h \
    bond.h \
    cell_list.h \
    cluster_analysis.h \
    data.h \
    force_field.h \
    rng.h \
    sim_box.h \
    types.h \
    xtcanalysis.h \
    welford.h \
    tensor_xyz.h \
    IO/io_input.h \
    IO/io_lammps.h \
    IO/io_pdb.h \
    xdrfile-1.1.4/include/xdrfile.h \
    xdrfile-1.1.4/include/xdrfile_trr.h \
    xdrfile-1.1.4/include/xdrfile_xtc.h \
    particles/particle_container.h \
    particles/chain.h \
    particles/cow.h \
    particles/dodecahedron.h \
    particles/ellipsoid.h \
    particles/lipid.h \
    particles/globular_sphere.h \
    particles/icosahedron.h \
    particles/oblatespheroid.h \
    particles/particle.h \
    particles/pentamer.h \
    particles/sphere.h \
    particles/spherepatch.h \
    particles/surface.h \
    particles/slab.h \
    particles/tennisball.h \
    Systems/flat_membrane.h \
    Systems/lipid_nanoparticle.h \
    Systems/membrane.h \
    Systems/system_base.h \
    Systems/system_container.h \
    Systems/vesicle.h \
    Systems/virus_pseudot3.h




