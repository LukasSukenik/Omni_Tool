#ifndef FLAT_MEMBRANE_H
#define FLAT_MEMBRANE_H

#include "system_base.h"
#include "atom.h"

#include "lipid.h"

#include "xtcanalysis.h"
#include "cluster_analysis.h"




class Cluster_Analysis : public System_Method
{
public:
    Cluster_Analysis() : System_Method() {}

    string help()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Flat_Membrane" << endl;
        ss << "System_execute: Cluster_Analysis" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "Atom_type: 2 3" << endl;
        ss << "Cluster_cutoff: 2.6" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate(Data& data)
    {
        data.in.param.validate_keyword("Input_type", "lammps_full");
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "traj_1.xtc");
        data.in.p_vec_int.validate_keyword("Atom_type", "1");
        data.in.p_float.validate_keyword("Cluster_cutoff", "2.6");
    }

    void execute(Data& data)
    {
        cerr << "Flat_Membrane::cluster_analysis" << endl;
        validate(data);

        Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        Clusters clusters(topo, data.in.p_vec_int["Atom_type"]); // list of particle indexes
        Trajectory traj(data.in.param["Trajectory_file"]);

        for(size_t i=0; i<traj.frame_count(); ++i)
        {
            topo.set_frame(traj[i]);
            cout << i << " " << clusters.analyze(topo, data.in.sim_box, data.in.p_float["Cluster_cutoff"]) << endl;
            clusters.clear();
        }
    }
};




class Detect_Pore : public System_Method
{
public:
    int number_of_cells;
    vector< vector<bool>> lattice;

    Detect_Pore() : System_Method() {}

    string help()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Flat_Membrane" << endl;
        ss << "System_execute: Calc_Pore" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "Atom_type: 2 3" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate(Data& data)
    {
        bool is_correct = true;
        data.in.p_float.validate_keyword("Cell_size", "0.15");

        if(data.in.p_vec_int["Atom_type"].empty())
        {
            cerr << "Atom_type empty" << endl;
            cerr << "add Atom_type: integer_1 interger_2 ... to input script, for example Atom_type: 2 3" << endl;
            is_correct = false;
        }

        if(!is_correct)
        {
            cerr << "Flat_Membrane::validate_pore_input exiting" << endl;
        }
    }

    void execute(Data& data)
    {
        cerr << "Flat_Membrane::is_Pore" << endl;
        validate(data);

        Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        Tensor_xyz box = data.in.sim_box.get_box();
        Trajectory traj;
        if(data.in.param.contains("Trajectory_file"))
        {
            traj.load(data.in.param["Trajectory_file"]);
        }



        cerr << "- cell size = " << data.in.p_float["Cell_size"] << endl;
        cerr << "- bead size = " << data.in.bead_size << endl;

        number_of_cells = box.x / data.in.p_float["Cell_size"];
        lattice = vector<vector<bool>>(number_of_cells, vector<bool>(number_of_cells, false));

        if(!data.in.param.contains("Trajectory_file"))
        {
            fill_lattice(data, topo, box);
            cout << is_pore() << endl; // true == 1
        }
        else
        {
            for(int i=0; i<traj.frame_count(); ++i)
            {
                topo.set_frame(traj[i]);
                box = traj.box_traj[i];

                fill_lattice(data, topo, box);
                cout << i << " " << is_pore() << endl; // true == 1
            }
        }
    }

    void set_lattice_false()
    {
        for(auto& row : lattice)
        {
            fill(row.begin(), row.end(), false);
        }
    }

    void fill_lattice(Data& data, Atoms& topo, Tensor_xyz& box)
    {
        double cell_size = box.x / number_of_cells;
        double inv_cell_size = number_of_cells / box.x; // box.x == box.y allways (enforced by lammps settings)

        bool do_iterations = (data.in.p_float["Cell_size"] < data.in.bead_size);
        int iterations = (0.5*data.in.bead_size / data.in.p_float["Cell_size"]) +1;
        double R2 = 0.5*data.in.bead_size * 0.5*data.in.bead_size;

        int cx=0, cy=0;
        int cell_ix=0, cell_iy=0;
        double cell_x_pos=0.0, cell_y_pos=0.0;
        double rx=0.0, ry=0.0;

        set_lattice_false(); // set all elements of latice to false -> false == no particle in cell

        for(Atom& a : topo)
        {
            for(int type : data.in.p_vec_int["Atom_type"])
            {
                if(a.type == type)
                {
                    //
                    // Object crosses cell border == cell is occupied
                    // - for every single cell ID "HIT" -> 4 cells are occupied [0,0][0,-1][-1,0][-1,-1] -> we are flooring the position to nearest cell
                    //
                    cx = binning_fce(a.pos.x, box.x, number_of_cells, inv_cell_size); // index of cell, where particle COM hit
                    cy = binning_fce(a.pos.y, box.x, number_of_cells, inv_cell_size);

                    if(do_iterations)
                    {
                        for(int dx=-iterations; dx<=iterations; ++dx) // iterate neighboring cells
                        {
                            for(int dy=-iterations; dy<=iterations; ++dy)
                            {
                                cell_ix = cx+dx;
                                cell_iy = cy+dy;
                                cell_x_pos = cell_ix*cell_size;
                                cell_y_pos = cell_iy*cell_size;
                                rx = cell_x_pos - (a.pos.x + 0.5*box.x); // distance of cell to particle + correct the position
                                ry = cell_y_pos - (a.pos.y + 0.5*box.y);

                                if(rx*rx + ry*ry <= R2) // only cells whose center fall inside the bead
                                {
                                    lattice[ wrap(cell_ix,   number_of_cells) ][ wrap(cell_iy,     number_of_cells) ] = true;
                                    lattice[ wrap(cell_ix,   number_of_cells) ][ wrap(cell_iy-1,   number_of_cells) ] = true;
                                    lattice[ wrap(cell_ix-1, number_of_cells) ][ wrap(cell_iy,     number_of_cells) ] = true;
                                    lattice[ wrap(cell_ix-1, number_of_cells) ][ wrap(cell_iy-1,   number_of_cells) ] = true;
                                }
                            }
                        }
                    }
                    else
                    {
                        lattice[cx][cy] = true;
                    }
                    break;
                }
            }
        }
    }

    bool is_pore()
    {
        for(int j=0; j<number_of_cells; ++j)
        {
            for(int k=0; k<number_of_cells; ++k)
            {
                if( lattice[j][k] == false )
                {
                    return true;
                }
            }
        }
        return false;
    }

    void debug(int i)
    {
        if(true || i==1107) // debug
        {
            for(int j=0; j<number_of_cells; ++j)
            {
                for(int k=0; k<number_of_cells; ++k)
                {
                    cout << ( lattice[j][k] ? "█" : " " );
                }
                cout << endl;
            }
            //std::cout << "\033[" << number_of_cells << "A";
        }
    }

    //
    // values are distributed from -box_size/2 to box_size/2
    //
    inline int binning_fce(double value, double box_size, int number_of_cells, double inv_cell_size)
    {
        if (value < -0.5*box_size) return 0;
        if (value >= 0.5*box_size) return number_of_cells - 1;

        return (value + 0.5*box_size) * inv_cell_size; // value / cell_size, cell_size = box/number_of_cells
    }

    //
    // x is from <0;count)
    // - with 0, without count
    //
    int wrap(int index, int count)
    {
        return (index + count) % count; // index can be -1, or >count
    }
};




class Flat_Membrane : public System_Base, public Particle
{
public:
    inline static const string keyword = "Flat_Membrane";
    const string name = "Flat_Membrane";

    Flat_Membrane() : System_Base("Flat_Membrane"), Particle("Flat_Membrane")  {}

    string help()
    {
        stringstream ss;

        ss << cluster_analysis.help() << endl;
        ss << detect_pore.help() << endl;

        return ss.str();
    }

    void execute(Data& data)
    {
        data.in.param.validate_keyword("System_execute", "Copy_Z | Calc_Z_Dist | Calc_Pore | Cluster_Analysis");
        if(data.in.param["System_execute"].compare("Copy_Z") == 0)           { copy_Z(data); }
        if(data.in.param["System_execute"].compare("Cluster_Analysis") == 0) { cluster_analysis.execute(data); }
        if(data.in.param["System_execute"].compare("Calc_Pore") == 0)        { detect_pore.execute(data); }
    }


    void generate( Data& data )
    {
        Lipids membrane = gen_flat_membrane(data);

        for(Lipid& lip : membrane)
        {
            beads.insert(beads.end(), lip.part.begin(), lip.part.end());
            bonds.insert(bonds.end(), lip.bond.begin(), lip.bond.end());
        }
    }

private:
    Cluster_Analysis cluster_analysis;
    Detect_Pore detect_pore;

    ///
    ///
    ///
    string help_gen_flat_membrane()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "Particle_type: Flat_Membrane" << endl;
        ss << "Num_lipids: 100" << endl;
        ss << "Number_of_receptors: 50" << endl;
        ss << "Mol_tag: 1" << endl;
        ss << "Leaflet_type: Upper Lower" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_gen_flat_membrane_input(Data& data)
    {
        data.in.p_int.validate_keyword("Num_lipids", "1000");
        data.in.p_int.validate_keyword("Number_of_receptors", "500");
        data.in.p_int.validate_keyword("Mol_tag", "1");
        data.in.p_int.validate_keyword("Receptor_type", "7");
        data.in.p_vec_string.validate_keyword("Leaflet_type", "Upper Lower");
        data.in.p_int.validate_keyword("ID", "1");
    }

    Lipids gen_flat_membrane(Data& data)
    {
        validate_gen_flat_membrane_input(data);

        int num_lipids = data.in.p_int["Num_lipids"];
        int num_receptors = data.in.p_int["Number_of_receptors"];
        int mol_tag = data.in.p_int["Mol_tag"];
        size_t size=0;
        if(!data.coll_beads.empty())
        {
            Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
            size = topo.size();
        }

        Lipids mem;
        int side_len = sqrt(num_lipids/2) +1;
        int count=0;
        double x=0.0,y=0.0;
        double z_up = 3.5, z_down=-3.5;
        double factor = pow(2.0, 1.0 / 6.0);

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

                    if( !exclude_XY(x,y,data) )
                    {
                        pos_up =   Tensor_xyz(x*factor, y*factor, z_up);
                        pos_down = Tensor_xyz(x*factor, y*factor, z_down);

                        mem.push_back(Lipid(pos_up,   dir_down, count,   mol_tag, get_leaflet_type( data.in.p_vec_string["Leaflet_type"][0] ) ));
                        mem.back().bond.offset(0, size);
                        mem.push_back(Lipid(pos_down, dir_up,   count+1, mol_tag, get_leaflet_type( data.in.p_vec_string["Leaflet_type"][1] ) ));
                        mem.back().bond.offset(0, size);
                        count+=2;
                    }
                }
            }
        }

        if(data.in.sim_box.xlo == 0.0)
        {
            data.in.sim_box.xlo = -0.5*side_len*factor + 0.5*factor;
            data.in.sim_box.xhi =  0.5*side_len*factor - 0.5*factor;
            data.in.sim_box.ylo = -0.5*side_len*factor + 0.5*factor;
            data.in.sim_box.yhi =  0.5*side_len*factor - 0.5*factor;
            data.in.sim_box.zlo = -0.5*side_len*factor + 0.5*factor;
            data.in.sim_box.zhi =  0.5*side_len*factor - 0.5*factor;
        }

        mem.set_receptor_type(data.in.p_int["Receptor_type"]);
        mem.convert_receptors(num_receptors);



        return mem;
    }

    Lipid::Leaflet get_leaflet_type(string type)
    {
        if(type.compare("Lower") == 0)
            return Lipid::Leaflet::lower;
        return Lipid::Leaflet::upper;
    }

    ///
    /// \brief exclude_XY - exclude point if its withing existing structure in XY plane
    /// \param data
    /// \return true if within structure
    ///
    bool exclude_XY(double x, double y, Data& data)
    {
        if(data.coll_beads.empty())
        {
            return false;
        }

        /*Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        double mil = 1000.0*1000.0;
        double patch = 0.5;
        Tensor_xyz o(x,y,0.0);

        for(Atom& a : topo)
        {
            if(a.pos.within(o, Tensor_xyz(patch, patch, mil)) ) // a is on the same X patch
            {
                return true;
            }
        }*/

        if(x*x + y*y < data.in.p_float["Exclude_radius"] * data.in.p_float["Exclude_radius"])
            return true;

        return false;
    }


    ///
    ///
    ///
    void copy_Z(Data& data)
    {
        cerr << "Flat_Membrane::execute -> Copy_Z" << endl;
        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& mem = data.coll_beads[sys_id];

        Lipids membrane_1 = Lipids(mem);
        Lipids membrane_2 = Lipids(mem, 2, mem.size()/4);

        membrane_1.move( Tensor_xyz(0.0, 0.0, data.in.system_var_a) );
        membrane_2.move( Tensor_xyz(0.0, 0.0, data.in.system_var_a*-1.0) );

        for(size_t i=0; i<mem.size()/4; ++i)
        {
            mem[4*i +0] = membrane_1[i].part[0];
            mem[4*i +1] = membrane_1[i].part[1];
            mem[4*i +2] = membrane_1[i].part[2];
            mem[4*i +3] = membrane_1[i].part[3];
        }

        for(Lipid& lip : membrane_2)
        {
            beads.insert(beads.end(), lip.part.begin(), lip.part.end());
            bonds.insert(bonds.end(), lip.bond.begin(), lip.bond.end());
        }

        data.coll_beads.push_back(beads);
        data.coll_bonds.push_back(bonds);

        cerr << "end of copy_Z" << endl;
    }

    int get_first_tail_bead(Atoms& a)
    {
        Lipid test;
        for(size_t i=0; i<a.size(); ++i)
        {
            if(test.is_tail(a[i].type))
            {
                return i;
            }
        }
        cerr << "Flat_Membrane::get_first_tail_bead - no bead identified as tail";
        exit(-1);
    }

};

#endif // FLAT_MEMBRANE_H
