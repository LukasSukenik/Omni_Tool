#ifndef LIPID_NANOPARTICLE_H
#define LIPID_NANOPARTICLE_H

#include "system_base.h"
#include "atom.h"
#include "xtcanalysis.h"

class Lipid_Nanoparticle : public System_Base
{
public:
    inline static const string keyword = "Lipid_Nanoparticle";
    const string name = "Lipid_Nanoparticle";

    Lipid_Nanoparticle() : System_Base("Lipid_Nanoparticle") {}

    void execute(Data& data)
    {
        if(data.in.system_function.compare("calc_water_content") == 0)
        {
            cerr << "calc_water_content::begin" << endl;

            // feeler atom
            Atom feeler;
            feeler.pos = Tensor_xyz(0.0, 0.0, 0.0);

            // box boundaries
            cerr << data.in.sim_box.xlo << " " << data.in.sim_box.xhi << endl;

            // trajectory load
            Trajectory traj;
            traj.load(data.in.trajectory);
            Atoms last_frame;
            last_frame.set_frame(traj[ traj.size()-1 ]);

            // feel stuff, similar to monte carlo pi calculation

            // example for a single location
            if( is_water(last_frame, feeler) )
            {
                cerr << "location is water" << endl;
            }
            else
            {
                cerr << "location is membrane" << endl;
            }

            // distribute location through box
            // count water and membrane location
            // profit?

            cerr << "calc_water_content::end" << endl;
        }
    }

    bool is_water(Atoms& frame, Atom& feeler)
    {
        return !frame.is_overlap(feeler, 1.0);
    }
};


#endif // LIPID_NANOPARTICLE_H
