#ifndef SIM_BOX_H
#define SIM_BOX_H

#include "types.h"
#include "atom.h"



class Simulation_Box
{
public:
    Simulation_Box() {}

    myFloat xlo = 0.0;
    myFloat xhi = 0.0;
    myFloat ylo = 0.0;
    myFloat yhi = 0.0;
    myFloat zlo = 0.0;
    myFloat zhi = 0.0;

    bool empty()
    {
        if( ( xlo == 0.0 && xhi == 0.0 ) || ( ylo == 0.0 && yhi == 0.0 ) || ( zlo == 0.0 && zhi == 0.0 ) )
                return true;
        return false;
    }

    inline double volume() {
        return (xhi-xlo)*(yhi-ylo)*(zhi-zlo);
    }

    Tensor_xyz get_random_pos()
    {
        Tensor_xyz a;
        a.x = ran() * (xhi-xlo) + xlo;
        a.y = ran() * (yhi-ylo) + ylo;
        a.z = ran() * (zhi-zlo) + zlo;
        return a;
    }

    bool is_in_box(Atoms atms)
    {
        for(auto& a : atms)
        {
            if(a.pos.x > xhi || a.pos.x < xlo)
                return false;
            if(a.pos.y > yhi || a.pos.y < ylo)
                return false;
            if(a.pos.z > zhi || a.pos.z < zlo)
                return false;
        }
        return true;
    }

    Atom usePBC(Tensor_xyz& pos_orig, double scale) const // used by dodecahedron only, dont touch
    {
        Tensor_xyz pos = pos_orig;

        while (pos.x < 0.0) {
            pos.x += xhi/scale;
        }
        while (pos.x > xhi/scale) {
            pos.x -= xhi/scale;
        }
        while (pos.y < 0.0) {
            pos.y += yhi/scale;
        }
        while (pos.y > yhi/scale) {
            pos.y -= yhi/scale;
        }
        while (pos.z < 0.0) {
            pos.z += zhi/scale;
        }
        while (pos.z > zhi/scale) {
            pos.z -= zhi/scale;
        }
        return pos;
    }

    friend std::istream& operator>>(std::istream& is, Simulation_Box& sim_box)
    {
        is >> sim_box.xlo >> sim_box.xhi >> sim_box.ylo >> sim_box.yhi >> sim_box.zlo >> sim_box.zhi;
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const Simulation_Box& box)
    {
        os << box.xlo << ", " << box.xhi << ", " << box.ylo << ", " << box.yhi << ", " << box.zlo << ", " << box.zhi;
        return os;
    }
};

#endif // SIM_BOX_H
