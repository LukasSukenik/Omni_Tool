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
        ss << help_traj_to_file() << endl;


        return ss.str();
    }

    void execute(Data& data)
    {
        data.in.param.validate_keyword("System_execute", "");

        // Fixing a gcmc xtc is not possible, xtc format cannot be written with variable per frame natoms
        //if(data.in.param["System_execute"].compare("Fix_GCMC_xtc") == 0) { fix_gcmc_xtc(data); }
        if(data.in.param["System_execute"].compare("Traj_to_file") == 0) { traj_to_file(data); }
    }

private:
    ///
    ///
    ///
    string help_file_to_traj()
    {
        stringstream ss;

        ss << "System_type: Generic" << endl;
        ss << "System_execute: File_to_traj" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_output_file: file.xtc" << endl;
        ss << "File_list: list" << endl;
        ss << "Only_last_frame:" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_file_to_traj_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_output_file", "file.xtc");
        data.in.param.validate_keyword("File_list", "list");
    }

    void file_to_traj(Data& data)
    {
        cerr << "Generic::file_to_traj" << endl;
        validate_file_to_traj_inputs(data);

        string in_file;
        Trajectory traj;
        IO_Lammps lammps;
        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];

        fstream fs(data.in.param["File_list"], fstream::in);
        if(!fs.is_open())
        {
            cout << "Failed to open file " << data.in.param["File_list"] << endl;;
        }
        else
        {
            while(!fs.eof())
            {
                fs >> in_file;

                lammps.load(in_file);
                traj.conf_traj.push_back(lammps.beads.get_all_pos());
                traj.box_traj.push_back( data.in.sim_box.get_box() );
                lammps.bonds.clear();
                lammps.angles.clear();
                lammps.beads.clear();
            }
            traj.write(data.in.param["Trajectory_output_file"]);
        }
        fs.close();
    }


    ///
    /// Traj last frame to file
    ///
    string help_traj_to_file()
    {
        stringstream ss;

        ss << "System_type: Generic" << endl;
        ss << "System_execute: Traj_to_file" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "Only_last_frame:" << endl;
        ss << "Output_type: pdb" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_traj_to_file_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "traj_1.xtc");
    }

    void traj_to_file(Data& data)
    {
        cerr << "Generic::traj_to_file" << endl;
        validate_traj_to_file_inputs(data);

        Trajectory traj(data);
        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];
        topo.set_frame(traj[0]);
    }

    ///
    /// Fix GCMC trajectory -> cant be done for xtc
    ///
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
