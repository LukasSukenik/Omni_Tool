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
    atom.h \
    chain.h \
    cow.h \
    data.h \
    dodecahedron.h \
    force_field.h \
    globular_sphere.h \
    icosahedron.h \
    input.h \
    oblatespheroid.h \
    particle.h \
    pentamer.h \
    rng.h \
    slab.h \
    sphere.h \
    spherepatch.h \
    surface.h \
    tennisball.h \
    xtcanalysis.h \
    welford.h \
    xdrfile-1.1.4/include/xdrfile.h \
    xdrfile-1.1.4/include/xdrfile_trr.h \
    xdrfile-1.1.4/include/xdrfile_xtc.h



