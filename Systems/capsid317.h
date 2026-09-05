#ifndef CAPSID317_H
#define CAPSID317_H

#include "system_base.h"
#include "atom.h"
#include "xtcanalysis.h"

class Capsid317 : public System_Base
{
public:
    inline static const string keyword = "Capsid317";
    const string name = "Capsid317";

    Capsid317() : System_Base("Capsid317") {}

    string help()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Capsid317" << endl;
        ss << "System_execute: Convert_Blender" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: md0001.xtc" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate(Data& data)
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "md0001.xtc");
        data.in.p_int.validate_keyword("Trajectory_frame", "1");
    }

    void execute(Data& data)
    {
        validate(data);
        if(data.in.param["System_execute"].compare("Convert_Blender") == 0) { Convert_Blender(data); }
    }

private:
    void Convert_Blender(Data& data)
    {
        /*cerr << "Capsid317::execute -> Convert_Blender" << endl;
        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& topo = data.coll_beads[sys_id];
        Trajectory traj(data);
        topo.set_frame(traj[ data.in.p_int["Trajectory_frame"] ]);

        Atoms penta[12];
        IO_XYZ penta_com;
        IO_XYZ penta_dir_up;
        IO_XYZ penta_dir_side;

        Atom temp1, temp2;

        for(int i=0; i<12; ++i)
        {
            penta_com.beads.push_back(topo.center_of_mass(i));

            penta[i] = topo.get_molecule(i);
            temp1 = penta[i].get_type(1).get_center_of_mass();
            temp2 = penta[i].get_type(5).get_center_of_mass();
            penta_dir_up.beads.push_back(temp1-temp2);

            temp1 = penta[i].center_of_mass(-1,0,36);
            temp2 = penta[i].center_of_mass(-1,36,71);
            penta_dir_side.beads.push_back(temp1-temp2);
        }

        penta_com.print_to_file("penta_coms");
        penta_dir_up.print_to_file("penta_dir_up");
        penta_dir_side.print_to_file("penta_dir_side");

        cerr << data.coll_beads.size() << " " << topo.size() << endl;
        //topo = topo.get_molecule(12); // remove pentamers, keep only chain
        cerr << data.coll_beads.size() << " " << topo.size() << endl;*/
    }
};

#endif // CAPSID317_H
