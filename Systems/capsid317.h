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
        cerr << "Capsid317::execute -> Convert_Blender" << endl;
        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& topo = data.coll_beads[sys_id];

        if( data.in.p_int["Trajectory_frame"] >= 0)
        {
            Trajectory traj(data);
            topo.set_frame(traj[ data.in.p_int["Trajectory_frame"] ]);
        }

        Atoms penta[12];
        IO_PDB penta_com;
        IO_PDB penta_normal;
        IO_PDB penta_tangens;

        Atom temp, normal, not_true_tangens, binormal, tangens, comm;
        size_t penta_size = 0;

        for(int i=0; i<12; ++i)
        {
            comm = topo.center_of_mass(i);

            comm.mol_tag = 1;

            penta_com.beads.push_back(comm);

            penta[i] = topo.get_molecule(i);
            temp = penta[i].get_type(1).get_center_of_mass();
            normal = temp-comm;
            normal.normalise();

            normal.mol_tag = 1;

            penta_normal.beads.push_back( normal );

            penta_size = penta[i].size();
            temp = penta[i].center_of_mass(i*penta_size + 0, i*penta_size + 37);
            not_true_tangens = temp-comm;
            not_true_tangens.normalise();
            binormal = normal.cross(not_true_tangens);
            binormal.normalise();
            tangens = binormal.cross(normal);
            tangens.normalise();

            tangens.mol_tag = 1;

            penta_tangens.beads.push_back(tangens);
        }

        penta_com.print_to_file("penta_coms.pdb");
        penta_normal.print_to_file("penta_normal.pdb");
        penta_tangens.print_to_file("penta_tangens.pdb");

        topo = topo.get_molecule(12); // remove pentamers, keep only chain
        topo.set_N(1);
    }
};

#endif // CAPSID317_H
