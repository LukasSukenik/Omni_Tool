#ifndef GENERIC_H
#define GENERIC_H

#include "system_base.h"
#include "xtcanalysis.h"




class Generic : public System_Base
{
public:
    inline static const string keyword = "Generic";
    const string name = "Generic";

    Generic() : System_Base("Generic") {}

    string help()
    {
        stringstream ss;

        //ss << help_fix_gcmc_xtc() << endl;

        return ss.str();
    }

    void execute(Data& data)
    {
        data.in.param.validate_keyword("System_execute", "");


        // Fixing a gcmc xtc is not possible, xtc format cannot be written with variable per frame natoms
        //if(data.in.param["System_execute"].compare("Fix_GCMC_xtc") == 0) { fix_gcmc_xtc(data); }
    }

private:
    string help_fix_gcmc_xtc()
    {
        stringstream ss;

        ss << "System_type: Generic" << endl;
        ss << "System_execute: Fix_GCMC_xtc" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.end" << endl;
        ss << "Trajectory_file: file.xtc" << endl;
        ss << "Trajectory_output_file: file2.xtc" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_fix_gcmc_xtc( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.end");
        data.in.param.validate_keyword("Trajectory_file", "file.xtc");
        data.in.param.validate_keyword("Trajectory_output_file", "file2.xtc");
    }

    // Fixing a gcmc xtc is not possible, xtc format cannot be written with variable per frame natoms
    void fix_gcmc_xtc(Data& data)
    {
        cerr << "Generic::fix_gcmc_xtc" << endl;
        validate_fix_gcmc_xtc(data);

        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];
        Trajectory traj;

        traj.load_gcmc_traj(data, topo);
        //traj.write(data.in.param["Trajectory_output_file"]);
    }

};

#endif // GENERIC_H
