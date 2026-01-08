#ifndef FLAT_MEMBRANE_H
#define FLAT_MEMBRANE_H

#include "system_base.h"
#include "atom.h"

#include "lipid.h"

class Flat_Membrane : public System_Base, public Particle
{
public:
    inline static const string keyword = "Flat_Membrane";
    const string name = "Flat_Membrane";

    Flat_Membrane() : System_Base("Flat_Membrane"), Particle("Flat_Membrane")  {}

    void execute(Data& data)
    {
        cerr << "flat membrane system" << endl;
    }

    void generate( Data& data )
    {
        Lipids membrane = gen_flat_membrane(data.in.num_lipids, data.in.num_rec, data.in.mol_tag);

        for(Lipid& lip : membrane)
        {
            beads.insert(beads.end(), lip.part.begin(), lip.part.end());
            bonds.insert(bonds.end(), lip.bond.begin(), lip.bond.end());
        }

        cerr << "suggesting box size: " << - 0.5*(sqrt(data.in.num_lipids/2)+1) << " to " << 0.5*(sqrt(data.in.num_lipids/2)+1) << endl;
    }

private:
    Lipids gen_flat_membrane(int num_lipids, int num_receptors, int mol_tag)
    {
        Lipids mem;

        int side_len = sqrt(num_lipids/2) +1;

        int count=0;
        double x=0.0,y=0.0;
        double z_up = 3.5, z_down=-3.5;

        Tensor_xyz pos_up = Tensor_xyz(0,0,z_up);
        Tensor_xyz pos_down = Tensor_xyz(0,0,z_down);

        Tensor_xyz dir_up = Tensor_xyz(0,0,1);
        Tensor_xyz dir_down = Tensor_xyz(0,0,-1);

        for(int i=0; i< side_len; ++i)
        {
            for(int j=0; j<side_len; ++j)
            {
                if(count < num_lipids)
                {
                    x = i - 0.5*side_len;
                    y = j - 0.5*side_len;

                    pos_up =   Tensor_xyz(x,y,z_up);
                    pos_down = Tensor_xyz(x,y,z_down);

                    mem.push_back(Lipid(pos_up,   dir_down, count,   mol_tag, Lipid::Leaflet::upper));
                    mem.push_back(Lipid(pos_down, dir_up,   count+1, mol_tag, Lipid::Leaflet::lower));
                    count+=2;
                }
            }
        }

        mem.convert_receptors(num_receptors);

        return mem;
    }

};

#endif // FLAT_MEMBRANE_H
