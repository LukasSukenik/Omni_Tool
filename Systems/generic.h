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
        ss << help_file_to_traj() << endl;

        return ss.str();
    }

    void execute(Data& data)
    {
        data.in.param.validate_keyword("System_execute", "");

        // Fixing a gcmc xtc is not possible, xtc format cannot be written with variable per frame natoms
        //if(data.in.param["System_execute"].compare("Fix_GCMC_xtc") == 0) { fix_gcmc_xtc(data); }
        if(data.in.param["System_execute"].compare("Traj_to_file") == 0) { traj_to_file(data); }
        if(data.in.param["System_execute"].compare("Coarse_grain") == 0) { coarse_grain(data); }
    }

private:
    ///
    ///
    ///
    string help_coarse_grain()
    {
        stringstream ss;

        ss << "System_type: Generic" << endl;
        ss << "System_execute: Coarse_grain" << endl;
        ss << "Input_type: pdb" << endl;
        ss << "Load_file: Single_Floor.pdb" << endl;
        ss << "Output_type: lammps_full" << endl;
        ss << "Center:" << endl; // centers the loaded pdb before running coarse_grain()
        ss << "Scale: 0.17" << endl;
        ss << "Atom_type: 8" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_coarse_grain_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.pdb");
    }

    void coarse_grain(Data& data)
    {
        cerr << "Generic::coarse_grain" << endl;
        validate_coarse_grain_inputs(data);

        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];
        Atoms temp;
        Tensor_xyz temp_pos;

        int i=1;
        for(int j=0; j<data.in.p_int["Number_of_floors"]; ++j)
        {
            for(Atom& a : topo)
            {
                if(a.atom_name.compare(" CA ") == 0)
                {
                    temp_pos = a.pos;
                    temp_pos.rotate(Tensor_xyz(0.0, 0.0, 1.0), j*deg_to_rad*data.in.p_float["z_rotation_deg"]);
                    temp.push_back( Atom(i, temp_pos + Tensor_xyz(0.0,0.0, j*data.in.p_float["z_dist"]), data.in.p_vec_int["Atom_type"][0], data.in.p_int["Mol_tag"]) );
                    ++i;
                }
            }
        }

        topo.clear();
        topo = temp;
        topo.move(data.in.p_tensor["Position_shift"]); // position shift is not applied for type system.
    }

    void calc_stuff() // only for bacteriophage project 18.5.2026
    {
        Tensor_xyz a_1(177.362000, 171.585007, 154.714996); // atoms from vmd to calculate floors z dist and rotation
        Tensor_xyz a_2(170.542007, 178.880997, 115.549004);
        Tensor_xyz com(193.50001525878906, 193.50001525878906, 125.25961303710938);

        double z_dist = a_1.z - a_2.z;
        a_1.z = 0.0;
        a_2.z = 0.0;
        com.z = 0.0;
        double angleDeg = (a_1-com).angleToDegrees(a_2-com);
        double angleRad = (a_1-com).angleTo(a_2-com);

        cerr << "z_dist: " << z_dist*0.13333333 << endl;
        cerr << "z_rotation_deg: " << angleDeg << endl;
    }

    ///
    /// Converts input files into snapshots of a trajectory
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
    /// Traj last frame of trajectory to standard output (file)
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
