#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H

#include "icosahedron.h"
#include "oblatespheroid.h"
#include "sphere.h"
#include "tennisball.h"
#include "spherepatch.h"
#include "pentamer.h"
#include "dodecahedron.h"
#include "chain.h"
#include "slab.h"
#include "globular_sphere.h"
#include "ellipsoid.h"
#include "vesicle.h"

using namespace std;

/**
 * @brief Structure_Container
 * map classes, each generating a different structure
 */
class Particle_Container : public map<string, Particle*>
{
public:
    Particle_Container()
    {
        //
        // Add structures
        //
        (*this)[Empty_Particle::keyword] = new Empty_Particle();
        (*this)[Icosahedron<Surface>::keyword] = new Icosahedron<Surface>();
        (*this)[Sphere::keyword] = new Sphere();
        (*this)[TennisBall::keyword] = new TennisBall();
        (*this)[OblateSpheroid::keyword] = new OblateSpheroid();

        (*this)[SpherePatch::keyword] = new SpherePatch();
        (*this)[Pentamer<PentaSurface>::keyword] = new Pentamer<PentaSurface>();
        (*this)[Dodecahedron::keyword] = new Dodecahedron();
        (*this)[Chain::keyword] = new Chain();
        (*this)[Slab::keyword] = new Slab();

        (*this)[NettedSlab::keyword] = new NettedSlab();
        (*this)[SphereJanus::keyword] = new SphereJanus();
        (*this)[Cow::keyword] = new Cow();
        (*this)[Globular_Sphere::keyword] = new Globular_Sphere();
        (*this)[Monomer::keyword] = new Monomer();

        (*this)[Ellipsoid::keyword] = new Ellipsoid();
        (*this)[Vesicle::keyword] = new Vesicle();
    }

    ~Particle_Container()
    {
        for (auto const& [key, val] : (*this))
        {
             delete val;
        }
    }

    void helpMessage()
    {
        cout << "Specify input files, run $ ./ico filename_1 filename_2 ...\n\n" << endl;
        cout << "Example input files:\n" << endl;

        for (auto& [key, val] : *this)
        {
            cout << (*val).help() << endl;
        }
        exit(1);
    }
};

#endif // PARTICLE_CONTAINER_H
