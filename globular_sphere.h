#ifndef MEHRNOOSH_SPHERE_H
#define MEHRNOOSH_SPHERE_H

#include "sphere.h"

class Globular_Sphere : public Sphere
{
public:
    inline static const string keyword = "globular_sphere";
    const string name = "globular_sphere";

    Globular_Sphere() : Sphere("Globular_Sphere") {}

    void generate( Data& data )
    {
        Atoms base_sphere;
        Atoms patch_sphere;
        Atoms full_sphere;

        // generate low count pseudoatom sphere, for repuls only interation
        fibonacci_sphere( base_sphere, data.in.num_of_beads, data.in.ff.lj[1].type, data.in.mol_tag); // size 0.4 beads

        // generate sphere for patches
        fibonacci_sphere( patch_sphere, data.in.num_lig, -1, data.in.mol_tag); // second fib. sphere, size 0.1 beads

        // scale patch_sphere by 1.0 + 1st sphere pseudoatom sigma - patch pseudoatom sigma
        patch_sphere.scale(1.0 + (data.in.ff.lj[1].sigma - data.in.ff.lj[2].sigma) / data.in.scale);

        // create the final particle based on defined patches
        for(auto& atom : patch_sphere)
        {
            for(auto& patch : data.in.patches)
            {
                if( (patch.x*(atom.x - patch.vx)) + (patch.y*(atom.y - patch.vy)) + (patch.z*(atom.z - patch.vz)) > 0 + (data.in.ff.lj[1].sigma - data.in.ff.lj[2].sigma) / data.in.scale )
                {
                    atom.type=patch.type;
                    full_sphere.push_back(atom);
                }
            }
        }

        bool add = true;
        for(auto& atom : base_sphere)
        {
            add = true;
            for(auto& patch : data.in.patches)
            {
                if( (patch.x*(atom.x - patch.vx)) + (patch.y*(atom.y - patch.vy)) + (patch.z*(atom.z - patch.vz)) > 0 )
                {
                    add = false;
                }
            }
            if(add)
            {
                full_sphere.push_back(atom);
            }
        }

        // insert the final structure to beads
        beads.insert(beads.end(), full_sphere.begin(), full_sphere.end());
    }
};

#endif // MEHRNOOSH_SPHERE_H
