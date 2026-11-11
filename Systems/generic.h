#ifndef GENERIC_H
#define GENERIC_H

#include "system_base.h"
#include "xtcanalysis.h"
#include "dsu.h"

#include <set>




class Make_Molecules_Whole : public System_Method
{
public:
    double h_b_x;
    double h_b_y;
    double h_b_z;
    Simulation_Box box;
    unordered_map<int, int> map_N; // [N] = index

    Make_Molecules_Whole() : System_Method() {}

    string help()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Generic" << endl;
        ss << "System_execute: Make_mol_whole" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate(Data& data)
    {
        data.in.param.validate_keyword("Input_type", "lammps_full");
        data.in.param.validate_keyword("Load_file", "data.start");
    }

    void execute(Data& data)
    {
        cerr << "Generic::make_molecules_whole" << endl;

        Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        Bonds& t_bo = data.coll_bonds[   data.id_map[ data.in.p_int["ID"] ]   ];
        box = data.in.sim_box;
        h_b_x = 0.5*(box.xhi - box.xlo);
        h_b_y = 0.5*(box.yhi - box.ylo);
        h_b_z = 0.5*(box.zhi - box.zlo);

        DSU dsu;
        unordered_map<int, set<int>> map_Bond; // bond [at1] = [at2], [at2 of second bond] && at1 < at2, ...

        for(Atom& a : topo) { box.set_central_image(a.pos); }
        map_N = topo.get_N_to_index_map();

        for(Bond& b : t_bo)
        {
            dsu.unite(b.at1, b.at2);
            map_Bond[ b.at1 ].insert( b.at2 );
            map_Bond[ b.at2 ].insert( b.at1 ); // every atom needs a list of all bonds
        }

        //for (const auto& [at1, at2_list] : map_Bond)
        for (const auto& [rN, molecule] : dsu.groups())
        {
            if(molecule.size() < 10)
            {
                for (int aN : molecule)
                {
                    make_whole(topo, rN, aN);
                }
            }
            else
            {
                dfs(topo, molecule, map_Bond);
            }
        }

        if(!is_molecules_whole(data))
        {
            cerr << "is_molecules_whole() failed" << endl;
            exit(-11);
        }
    }

    void make_rigid_whole()
    {
        /*cerr << "Generic::make_rigid_molecules_whole" << endl;
        cerr << "Asuming phage system, with cell wall types 8 9" << endl;
        cerr << "Works by assuming neighbors in incremental index" << endl;

        //Tensor_xyz pbc, pbc_inv;
        //pbc.x = box.xhi - box.xlo;
        //pbc.x = box.yhi - box.ylo;
        //pbc.x = box.zhi - box.zlo;

        //pbc_inv.x = 1.0 / pbc.x;
        //pbc_inv.y = 1.0 / pbc.y;
        //pbc_inv.z = 1.0 / pbc.z;
        //a.distSQ_pbc(a, pbc, pbc_inv);

        //Info) index: 45657
        //Info) x: 8.306784
        //Info) y: 15.927391
        //Info) z: 20.920244
        size_t origin = 45657; // topo index, N=45658 -> should be the same for all

        size_t curr = origin;
        for(size_t i=origin-1; i>=0; ++i)
        {
            if(topo[i].type == 8 || topo[i].type == 9)
            {
                make_whole(topo, topo[curr].N, topo[i].N);
                curr = i;
            }
        }

        curr = origin;
        for(size_t i=origin+1; i<topo.size(); ++i)
        {
            if(topo[i].type == 8 || topo[i].type == 9)
            {
                make_whole(topo, topo[curr].N, topo[i].N);
                curr = i;
            }
        }*/
    }

    void make_whole(Atoms& topo, int rN, int aN)
    {
        int i;
        int j;
        Tensor_xyz dir;
        Tensor_xyz fix;
        Tensor_xyz o;

        if(rN != aN)
        {
            i = map_N[ rN ];
            j = map_N[ aN ];

            o = topo[j].pos;

            dir = topo[i].pos - topo[j].pos;
            fix.x = (abs(dir.x) > h_b_x) ? box.xhi-box.xlo : 0.0;
            fix.y = (abs(dir.y) > h_b_y) ? box.yhi-box.ylo : 0.0;
            fix.z = (abs(dir.z) > h_b_z) ? box.zhi-box.zlo : 0.0;

            //  dir > 0.0 ? (+i;-j) : (+i;+j)    dir < 0.0 ? (-i;+j) : (-i;-j)
            topo[j].pos.x += (dir.x > 0.0) ? fix.x : -fix.x;
            topo[j].pos.y += (dir.y > 0.0) ? fix.y : -fix.y;
            topo[j].pos.z += (dir.z > 0.0) ? fix.z : -fix.z;

            dir = topo[i].pos - topo[j].pos;
            if(abs(dir.x) > h_b_x || abs(dir.y) > h_b_y || abs(dir.z) > h_b_z)
            {
                cerr << topo[i].pos << endl;
                cerr << o << endl;
                cerr << "->" << endl;
                cerr << topo[j].pos << endl;
                exit(1);
            }
        }
    }

    ///
    /// Depth First graph algo
    ///
    void dfs(Atoms& topo,
             const set<int>& molecule,
             unordered_map<int, set<int>>& map_Bond)
    {
        if(molecule.empty())
            return;

        unordered_set<int> seen;
        vector<int> stack;

        int curr = *molecule.begin();
        stack.push_back(curr);

        while(!stack.empty())
        {
            curr = stack.back();
            stack.pop_back();

            if(seen.contains(curr))
                continue;      // already visited, skip
            seen.insert(curr);

            if(map_Bond.contains(curr)) // true if to found
            {
                for(int next : map_Bond[curr] ) // loop over set<int>>&
                {
                    if(!seen.contains(next))
                    {
                        make_whole(topo, curr, next);
                        stack.push_back(next);
                    }
                }
            }
        }
        if(molecule.size() != seen.size())
        {
            cerr << "Molecules size: " << molecule.size() << ", seen: " << seen.size() << endl;
            exit(-1);
        }
    }

    bool is_molecules_whole(Data& data)
    {
        Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        Bonds& t_bo = data.coll_bonds[   data.id_map[ data.in.p_int["ID"] ]   ];
        int i,j;

        for(Bond& b : t_bo)
        {
            i = map_N[b.at1];
            j = map_N[b.at2];

            if(topo[i].dist(topo[j]) > h_b_x)
            {
                cerr << "Bond huge " << b.at1 << " : " << b.at2 << " = " << topo[i].dist(topo[j]) << endl;
                return false;
            }
        }
        return true;
    }
};




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
        ss << help_load_traj() << endl;
        ss << help_file_to_traj() << endl;

        return ss.str();
    }

    void execute(Data& data)
    {
        data.in.param.validate_keyword("System_execute", "");

        // Fixing a gcmc xtc is not possible, xtc format cannot be written with variable per frame natoms
        //if(data.in.param["System_execute"].compare("Fix_GCMC_xtc") == 0) { fix_gcmc_xtc(data); }
        if(data.in.param["System_execute"].compare("Load_traj") == 0) { load_traj(data); }
        if(data.in.param["System_execute"].compare("Coarse_grain") == 0) { coarse_grain(data); }
        if(data.in.param["System_execute"].compare("Get_boundaries") == 0) { get_boundaries(data); }
        if(data.in.param["System_execute"].compare("Make_mol_whole") == 0) { m_M_W.execute(data); }
    }

private:
    Make_Molecules_Whole m_M_W;

    string help_get_boundaries()
    {

    }

    void validate_get_boundaries_inputs( Data& data )
    {
        if(data.in.p_vec_int["Atom_type"].empty())
        {
            cerr << "Atom_type empty" << endl;
            cerr << "add Atom_type: integer_1 interger_2 ... to input script, for example Atom_type: 2 3" << endl;
            exit(-1);
        }
    }

    void get_boundaries(Data& data)
    {
        cerr << "Generic::get_boundaries" << endl;

        Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        Atoms min20, max20;
        min20.reserve(20);
        max20.reserve(20);

        for(Atom& a : topo)
        {
            if(a.type == data.in.p_vec_int["Atom_type"][0])
            {
                if(min20.size() < 20) { min20.push_back(a); }
                else
                    replace_max_w_min(min20, a);

                if(max20.size() < 20) { max20.push_back(a); }
                else
                    replace_min_w_max(max20, a);
            }
        }

        cout << average_z(min20) << " " << average_z(max20) << endl;
    }

    void replace_max_w_min(Atoms& o, Atom a)
    {
        double min=999.9;
        int index=-1;

        for(int i=0; i<o.size(); ++i)
        {
            if(min > o[i].pos.z && o[i].pos.z > a.pos.z)
            {
                min = o[i].pos.z;
                index = i;
            }
        }

        if(index != -1)
        {
            o[index].pos.z = a.pos.z;
        }
    }

    void replace_min_w_max(Atoms& o, Atom a)
    {
        double max=-999.9;
        int index=-1;

        for(int i=0; i<o.size(); ++i)
        {
            if(max < o[i].pos.z && o[i].pos.z < a.pos.z)
            {
                max = o[i].pos.z;
                index = i;
            }
        }

        if(index != -1)
        {
            o[index].pos.z = a.pos.z;
        }
    }

    double average_z(Atoms& o)
    {
        double sum_z = 0.0;
        for(Atom& a : o)
        {
            sum_z += a.pos.z;
        }
        return sum_z / o.size();
    }

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
                if(strcmp(a.atom_name, " CA ") == 0)
                {
                    temp_pos = a.pos;
                    temp_pos.rotate(Tensor_xyz(0.0, 0.0, 1.0), j*deg_to_rad*data.in.p_float["z_rotation_deg"]);
                    temp.push_back( Atom(i, temp_pos + Tensor_xyz(0.0,0.0, j*data.in.p_float["z_dist"]), data.in.p_vec_int["Atom_type"][0], data.in.p_int["Mol_tag"]) );
                    ++i;
                }
            }
        }

        double radius = 4.0;
        double z_pos = 0.0;// -0.5* data.in.p_float["z_dist"] * data.in.p_int["Number_of_floors"];
        for(int x=-5; x<6; ++x)
        {
            for(int y=-5; y<6; ++y)
            {
                if(x*x + y*y < radius*radius)
                {
                    temp_pos = Tensor_xyz(x,y,z_pos);
                    temp.push_back( Atom(i, temp_pos, data.in.p_vec_int["Atom_type"][0], data.in.p_int["Mol_tag"]) );
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




    string help_load_traj()
    {
        stringstream ss;

        ss << "System_type: Generic" << endl;
        ss << "System_execute: Load_traj" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "Trajectory_frame: 5" << endl;
        ss << "Output_type: pdb" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_load_traj_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "traj_1.xtc");
        data.in.p_int.validate_keyword("Trajectory_frame", "5");
    }

    void load_traj(Data& data)
    {
        validate_load_traj_inputs(data);

        Trajectory traj(data);
        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];
        topo.set_frame(traj[ data.in.p_int["Trajectory_frame"] ]);
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
