#ifndef LIPID_H
#define LIPID_H

#include <vector>

#include "atom.h"
#include "bond.h"



#define RECEPTOR_BOND_TYPE 4
#define RECEPTOR_ANGLE_TYPE 1

#define TAIL_HEAD_BOND 1
#define TAIL_TAIL_BOND 2
#define HARMONIC_BOND 3
#define TAIL_END_BOND 4
#define HARMONIC_BOND2 5

using namespace std;



class Lipid {
public:
    Lipid(){}

    enum class Leaflet {upper, lower};

    static const int head_upper_leaf = 1;
    static const int tail_upper_leaf = 2;
    static const int tail_2_upper_leaf = 3;
    static const int head_lower_leaf = 4;
    static const int tail_lower_leaf = 5;
    static const int tail_2_lower_leaf = 6;
    static const int receptor = 7;

    Atoms part;  // consecutive:  0:HEAD, 1:TAIL, 2:TAIL (the further one)
    Bonds bond;      // 0: between 1-2 TAIL_HEAD_BOND, 1: 2-3 TAIL_TAIL_BOND, 2: 1-3 HARMONIC_BOND, 3: 3-4 TAIL_END_BOND, 4: 2-4 HARMONIC_BOND2

    Lipid(Atom p1, Atom p2, Atom p3, Atom p4)
    {
        part.resize(4);
        bond.resize(5);

        part[0] = p1;
        part[1] = p2;
        part[2] = p3;
        part[3] = p4;

        changeN_part(p1.N, 1);
        changeN_bond(p1.N);

        bond[0].type = TAIL_HEAD_BOND;
        bond[1].type = TAIL_TAIL_BOND;
        bond[2].type = HARMONIC_BOND;
        bond[3].type = TAIL_END_BOND;
        bond[4].type = HARMONIC_BOND2;
    }

    Atom get_direction()
    {
        Atom dir = part[0] - part[3]; // head - tail
        dir.normalise();
        return dir;
    }

    void changeN_part(int N, int mol_tag) {
        part[0].N = N;
        part[1].N = N+1;
        part[2].N = N+2;
        part[3].N = N+3;

        part[0].mol_tag = mol_tag;
        part[1].mol_tag = mol_tag;
        part[2].mol_tag = mol_tag;
        part[3].mol_tag = mol_tag;

    }

    void changeN_bond(int N) {
        bond[0].N = N;
        bond[1].N = N+1;
        bond[2].N = N+2;
        bond[3].N = N+3;
        bond[4].N = N+4;

        bond[0].at1 = N;
        bond[0].at2 = N+1;

        bond[1].at1 = N+1;
        bond[1].at2 = N+2;

        bond[2].at1 = N;
        bond[2].at2 = N+2;

        bond[3].at1 = N+2;
        bond[3].at2 = N+3;

        bond[4].at1 = N+1;
        bond[4].at2 = N+3;
    }

    void set_bead_type(Leaflet leaf)
    {
        if(leaf == Leaflet::upper)
        {
            part[0].type = head_upper_leaf;
            part[1].type = tail_upper_leaf;
            part[2].type = tail_upper_leaf;
            part[3].type = tail_2_upper_leaf;
        }
        if(leaf == Leaflet::lower)
        {
            part[0].type = head_lower_leaf;
            part[1].type = tail_lower_leaf;
            part[2].type = tail_lower_leaf;
            part[3].type = tail_2_lower_leaf;
        }
    }
};




class Lipids : public vector<Lipid>
{
public:
    bool update_positions(vector<Tensor_xyz>& frame)
    {
        if(4*size() != frame.size())
        {
            return false;
        }

        for(unsigned int i=0; i<size(); ++i) // loop over lipids
        {
            Lipid& lip = this->at(i);
            lip.part[0].pos = frame[4*i +0];
            lip.part[1].pos = frame[4*i +1];
            lip.part[2].pos = frame[4*i +2];
            lip.part[3].pos = frame[4*i +3];
        }

        return true;
    }

    Atoms get_Atoms()
    {
        Atoms at;
        for(Lipid& lip : (*this))
        {
            at.insert(at.end(), lip.part.begin(), lip.part.end());
        }
        return at;
    }
};




bool isTop(int i, vector<Atom>& particles) {
    //
    // search particles, if another particle of type HEAD found under i, then i is top layer
    //
    for(unsigned int j=0; j<particles.size(); ++j) {
        if(particles[j].type == Lipid::head_upper_leaf ||  particles[j].type == Lipid::head_lower_leaf) {
            if( ( (particles[j].pos.x - particles[i].pos.x)*(particles[j].pos.x - particles[i].pos.x)
                  + (particles[j].pos.y - particles[i].pos.y)*(particles[j].pos.y - particles[i].pos.y)
                  + (particles[j].pos.z - particles[i].pos.z)*(particles[j].pos.z - particles[i].pos.z) ) < 36.0 ) {
                if( (particles[j].pos.z - 4.0 > particles[i].pos.z) ) {
                    return true;
                }
            }
        }
    }
    return false;
}

#endif // LIPID_H
