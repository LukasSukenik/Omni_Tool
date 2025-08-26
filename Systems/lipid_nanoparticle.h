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

    Force_Field ff;

    Lipid_Nanoparticle() : System_Base("Lipid_Nanoparticle") {}

    void execute(Data& data)
    {
        if(data.in.system_function.compare("calc_water_content") == 0)
        {
            cerr << "calc_water_content::begin" << endl;

            // feeler atom
            Atom feeler;
            feeler.pos = Tensor_xyz(0.0, 0.0, 0.0);
            feeler.type = 2; // lipid tail

            // set ff
            ff.lj[1] = LJ(-0.8);
            ff.lj[2] = LJ(2.8);
            ff.lj[3] = LJ(2.8);

            // box boundaries
            cerr << data.in.sim_box.xlo << " " << data.in.sim_box.xhi << endl;

            // trajectory load
            Trajectory traj;
            traj.load(data.in.trajectory);

            int sys_id = data.id_map[data.in.id];
            Atoms& last_frame = data.coll_beads[sys_id];
            last_frame.set_frame(traj[ traj.size()-1 ]);
            cerr << "Trajectory loaded" << endl;

            double step_size = 1.0;
            int water_count = 0;
            int membrane_count = 0;

            double range = data.in.sim_box.xhi - data.in.sim_box.xlo;
            for(double x=data.in.sim_box.xlo; x<data.in.sim_box.xhi; x+= step_size) {
                cerr << "\r" << 100.0 * (x - data.in.sim_box.xlo) / range << " %         " << flush;
                for (double y=data.in.sim_box.ylo; y<data.in.sim_box.yhi; y+=step_size) {
                    for (double z=data.in.sim_box.zlo; z<data.in.sim_box.zhi; z+= step_size) {
                        feeler.pos = Tensor_xyz(x, y, z);
                        if( is_water(last_frame, feeler) )
                        {
                            water_count++;
                        }
                        else
                        {
                            membrane_count++;
                        }
                    }
                }
            }

            if (membrane_count != 0) {
                double ratio = static_cast<double>(water_count) / membrane_count;
                cout << "Water to Membrane ratio:" << ratio << endl;
                cout << "Membrane %:" << 100.0 * membrane_count / (membrane_count + water_count) << endl;
                cout << "step size:" << step_size << endl;
            } else {
                cerr << "Membrane count is zero, cannot compute ratio." << endl;
            }

            cerr << "calc_water_content::end" << endl;
        }
    }

    bool is_water(Atoms& frame, Atom& feeler)
    {
        return !frame.is_overlap(feeler, 1.0);
    }
};


#endif // LIPID_NANOPARTICLE_H
