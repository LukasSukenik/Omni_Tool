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
    }

    void execute(Data& data)
    {
        validate(data);
        if(data.in.param["System_execute"].compare("Convert_Blender") == 0) { Convert_Blender(data); }
    }

private:
    IO_PDB rna_com = IO_PDB("rna_coms.pdb");
    IO_PDB rna_normal = IO_PDB("rna_normals.pdb");
    IO_PDB rna_tangens = IO_PDB("rna_tangens.pdb");

    IO_PDB penta_com = IO_PDB("penta_coms.pdb");
    IO_PDB penta_normal = IO_PDB("penta_normals.pdb");
    IO_PDB penta_tangens = IO_PDB("penta_tangens.pdb");

    Trajectory traj_rna_com = Trajectory("", "rna_coms.xtc");
    Trajectory traj_rna_normal = Trajectory("", "rna_normals.xtc");
    Trajectory traj_rna_tangens = Trajectory("", "rna_tangens.xtc");

    Trajectory traj_penta_com = Trajectory("", "penta_coms.xtc");
    Trajectory traj_penta_normal = Trajectory("", "penta_normals.xtc");
    Trajectory traj_penta_tangens = Trajectory("", "penta_tangens.xtc");




    void Convert_Blender(Data& data)
    {
        cerr << "Capsid317::execute -> Convert_Blender" << endl;
        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& topo = data.coll_beads[sys_id];
        Trajectory traj(data);

        if( data.in.p_vec_int.contains("Trajectory_settings") )
        {
            for(size_t i=0; i<traj.frame_count(); ++i)
            {
                topo.set_frame( traj[i] );
                build_single_frame(topo);
                if(i==0) { write_pdbs(); }
                if(i>0)  { store_xtcs(i, topo, traj); }
            }
            write_xtcs(traj.get_precision());
        }

        if( data.in.p_int.contains("Trajectory_frame") )
        {
            topo.set_frame(traj[ data.in.p_int["Trajectory_frame"] ]);
            build_single_frame(topo);
            write_pdbs();
        }
    }

    void build_single_frame(Atoms& topo)
    {
        penta_com.beads = build_pentas_com(topo);
        penta_normal.beads = build_pentas_normals(topo);
        penta_tangens.beads = build_pentas_tangens(topo);

        rna_com.beads = build_RNA_chain(topo);
    }

    void write_pdbs()
    {
        penta_com.print_to_file();
        penta_normal.print_to_file();
        penta_tangens.print_to_file();

        rna_com.print_to_file();
    }

    void write_xtcs(float precision)
    {
        traj_rna_com.set_precision( precision );
        traj_rna_com.write();
    }

    void store_xtcs(size_t i, Atoms& topo, Trajectory& traj)
    {
        traj_penta_com.store_single_frame(i, traj.time_traj[i], traj.box_traj[i], penta_com.beads.get_all_pos());
        traj_penta_normal.store_single_frame(i, traj.time_traj[i], traj.box_traj[i], penta_normal.beads.get_all_pos());
        traj_penta_tangens.store_single_frame(i, traj.time_traj[i], traj.box_traj[i], penta_tangens.beads.get_all_pos());

        traj_rna_com.store_single_frame(i, traj.time_traj[i], traj.box_traj[i], rna_com.beads.get_all_pos());
    }

    Atoms build_RNA_chain(Atoms& topo)
    {
        Atoms rna = topo.get_molecule(12); // remove pentamers, keep only chain
        rna.set_N(1);
        return rna;
    }

    Atoms build_pentas_com(Atoms& topo)
    {
        Atom comm;
        Atoms p_coms;
        for(int i=0; i<12; ++i)
        {
            comm = topo.center_of_mass(i);
            comm.mol_tag = 1;
            p_coms.push_back(comm);
        }
        return p_coms;
    }

    Atoms build_pentas_normals(Atoms& topo)
    {
        Atom temp, normal, comm;
        Atoms p_normals;
        for(int i=0; i<12; ++i)
        {
            comm = topo.center_of_mass(i);

            temp = topo.get_molecule(i).get_type(1).get_center_of_mass();
            normal = temp-comm;
            normal.normalise();
            normal.mol_tag = 1;

            p_normals.push_back( normal );
        }
        return p_normals;
    }

    Atoms build_pentas_tangens(Atoms& topo)
    {
        Atoms penta[12];
        Atom temp, normal, not_true_tangens, binormal, tangens, comm;
        size_t penta_size = 0;
        Atoms p_tangens;
        for(int i=0; i<12; ++i)
        {
            comm = topo.center_of_mass(i);

            penta[i] = topo.get_molecule(i);
            temp = penta[i].get_type(1).get_center_of_mass();
            normal = temp-comm;
            normal.normalise();
            normal.mol_tag = 1;

            penta_size = penta[i].size();
            temp = penta[i].center_of_mass(i*penta_size + 0, i*penta_size + 37);
            not_true_tangens = temp-comm;
            not_true_tangens.normalise();
            binormal = normal.cross(not_true_tangens);
            binormal.normalise();
            tangens = binormal.cross(normal);
            tangens.normalise();
            tangens.mol_tag = 1;

            p_tangens.push_back(tangens);
        }
        return p_tangens;;
    }
};

#endif // CAPSID317_H
