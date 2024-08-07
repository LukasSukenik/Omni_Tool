TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += xdrfile-1.1.4/include/

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
    chain.h \
    cow.h \
    data.h \
    dodecahedron.h \
    ellipsoid.h \
    force_field.h \
    globular_sphere.h \
    icosahedron.h \
    io_input.h \
    io_lammps.h \
    io_pdb.h \
    oblatespheroid.h \
    particle.h \
    pentamer.h \
    rng.h \
    sim_box.h \
    slab.h \
    sphere.h \
    spherepatch.h \
    surface.h \
    tennisball.h \
    types.h \
    xtcanalysis.h \
    welford.h \
    xdrfile-1.1.4/include/xdrfile.h \
    xdrfile-1.1.4/include/xdrfile_trr.h \
    xdrfile-1.1.4/include/xdrfile_xtc.h



