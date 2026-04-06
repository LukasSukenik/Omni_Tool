#ifndef FLAT_MEMBRANE_H
#define FLAT_MEMBRANE_H

#include "system_base.h"
#include "atom.h"

#include "lipid.h"

#include "xtcanalysis.h"
#include "cluster_analysis.h"






class Flat_Membrane : public System_Base, public Particle
{
public:
    inline static const string keyword = "Flat_Membrane";
    const string name = "Flat_Membrane";

    Flat_Membrane() : System_Base("Flat_Membrane"), Particle("Flat_Membrane")  {}

    string help()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "Particle_type: Flat_Membrane\n";
        ss << "System_type: Flat_Membrane" << endl;

        return ss.str();
    }

    void execute(Data& data)
    {
        cerr << "Flat_Membrane::execute" << endl;
        validate_load_file(data);
        validate_ID(data);

        int sys_id = data.id_map[ data.in.param_int["ID"] ];
        Atoms& mem = data.coll_beads[sys_id];
        vector<int> mol_tags = mem.get_Mol_Types();

        if(data.in.system_function.compare("Copy_Z") == 0)           { copy_Z(data, mem); }
        if(data.in.system_function.compare("Cluster_Analysis") == 0) { cluster_analysis(data, mem); }
        if(data.in.system_function.compare("Calc_Z_Dist") == 0)      { calc_Z_Dist(data, mem); }
        if(data.in.system_function.compare("Calc_Pore") == 0)        { is_Pore(data, mem); }
    }


    void generate( Data& data )
    {
        Lipids membrane = gen_flat_membrane(data.in.param_int["Num_lipids"], data.in.param_int["Number_of_receptors"], data.in.param_int["Mol_tag"]);

        for(Lipid& lip : membrane)
        {
            beads.insert(beads.end(), lip.part.begin(), lip.part.end());
            bonds.insert(bonds.end(), lip.bond.begin(), lip.bond.end());
        }

        cerr << "suggesting box size: " << - 0.5*(sqrt(data.in.param_int["Num_lipids"]/2)+1) << " to " << 0.5*(sqrt(data.in.param_int["Num_lipids"]/2)+1) << endl;
    }

private:
    void validate_select_types(Data& data)
    {
        if(!data.in.param_vector_int.contains("Atom_type"))
        {
            cerr << "Missing keyword; Atom_type: 2 3" << endl;
            exit(-1);
        }
        if(data.in.param_vector_int["Atom_type"].empty())
        {
            cerr << "Atom_type: empty; Atom_type: 2 3" << endl;
            exit(-1);
        }
    }

    void validate_cluster_cutoff(Data& data)
    {
        if(!data.in.param_float.contains("Cluster_cutoff"))
        {
            cerr << "Missing keyword; Cluster_cutoff: 2.6" << endl;
            exit(-1);
        }
        if(data.in.param_float["Cluster_cutoff"] < 0.01 && data.in.param_float["Cluster_cutoff"] > 5.0)
        {
            cerr << "WARNING: Cluster_cutoff: is set to " << data.in.param_float["Cluster_cutoff"] << endl;
        }
    }

    void validate_load_file(Data& data)
    {
        if(!data.in.param.contains("Load_file"))
        {
            cerr << "Missing keyword; Load_file: data.start" << endl;
            exit(-1);
        }
    }

    void validate_ID(Data& data)
    {
        if(!data.in.param_int.contains("ID"))
        {
            cerr << "Missing keyword; Load_file: data.start" << endl;
            cerr << "potentially Missing keyword; ID: 0, set to 0 by default during Data::load_data()" << endl;
            exit(-1);
        }
    }

    void validate_trajectory(Data& data)
    {
        if(!data.in.param.contains("Trajectory_file"))
        {
            cerr << "Missing keyword; Trajectory_file: traj_1.xtc" << endl;
            exit(-1);
        }
    }




    void calc_Z_Dist(Data& data, Atoms& mem)
    {
        cerr << endl;
        cerr << "Calc Z Dist:" << endl;

        //
        // Calculate distance between 2 peaks of histogram in Z axis
        //
        // Histogram is an array of integers initialized to 0, use vector class
        // - https://www.w3schools.com/cpp/cpp_arrays.asp
        // - https://www.w3schools.com/cpp/cpp_for_loop.asp
        // - https://www.w3schools.com/cpp/cpp_operators.asp
        // - https://www.w3schools.com/cpp/cpp_vectors.asp
        //
        // The for loop over particles to determine their Z coordinate and increment the corresponding element in the histogram
        //
        // Based on the histogram, determine the distance of the 2 membranes - whether they have stalk or not
        // - dont bother with periodic boundary conditions at first, use the pulling simulations which do not have moment across the periodic boundary
        //
        // Take advantage of LLM, chatGPT, Claude are great for this
        // - but the most important skill here is learning algorithmic thinking and procedural decomposition
        // -- i.e. thinking like a computer
        // --- This will realy help later down the line when you will be analyzing your simulations
        // ---- you don't want to be the guy who manually types in indexes for energy calculations in gromacs for every simulation
        //
        vector<int> histogram;
        Trajectory traj;
        traj.load(data.in.param["Trajectory_file"]);

        for(int i=0; i<traj.frame_count(); ++i) // looping over trajectory frames
        {
            mem.set_frame(traj[i]);
        }
        // or you can access the trajectory directly
        cerr << "frame 0, particle 1: " << traj[0][1].x << endl;

        int atom_id = 0;
        cerr << "accessing atom z coordinate - pos = position: " << mem[atom_id].pos.z << endl;
        cerr << "Total atom count: " << mem.size() << " or " << traj[atom_id].size() << endl;

        //
        // You need a binning function for that
        // - convert continuous Z coordinate of floating point format to discrete integers value for histogram element index
        // - look up floor function,
        //
        int a = (int) floor(0.57);
        cerr << "a: " << a << endl;

        //
        // Now we just determine the distance based on the histogram
        //

        cerr << endl;
    }




    /**
     * @brief cluster_analysis
     * @param data
     * @param mem
     *
     * System_type: Flat_Membrane       -- needed call method cluster_analysis
     * System_execute: Cluster_Analysis -- needed call method cluster_analysis
     * Select_types: 1             -> validate
     * Cluster_cutoff: 2.6         -> validate
     * Input_type: lammps_full     -- input type does not matter, except none
     * Load_file: data.start       -- validated in Flat_Membrane::execute
     * Trajectory_file: traj_1.xtc -> validate
     * ID: 1                       -- validated in Flat_Membrane::execute
     *
     */
    void cluster_analysis(Data& data, Atoms& mem)
    {
        cerr << "Flat_Membrane::cluster_analysis" << endl;

        validate_select_types(data);
        validate_cluster_cutoff(data);
        validate_trajectory(data);

        for(int type : data.in.param_vector_int["Atom_type"])
        {
            cerr << type << ", ";
        }
        cerr << " | cutoff: " << data.in.param_float["Cluster_cutoff"] << endl;;

        Clusters clusters(mem, data.in.param_vector_int["Atom_type"]); // list of particle indexes
        Trajectory traj;
        traj.load(data.in.param["Trajectory_file"]);

        for(int i=0; i<traj.frame_count(); ++i)
        {
            mem.set_frame(traj[i]);
            cout << i << " " << clusters.analyze(mem, data.in.param_float["Cluster_cutoff"]) << endl;
            clusters.clear();
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

    bool validate_pore_input(Data& data)
    {
        bool is_correct = true;
        if(data.in.cell_size == 0.0)
        {
            cerr << "Cell_size set to 0.0 by default" << endl;
            cerr << "add Cell_size: floating_point_value to input script, for example 0.15" << endl;
            is_correct = false;
        }

        if(data.in.param_vector_int["Atom_type"].empty())
        {
            cerr << "Atom_type empty" << endl;
            cerr << "add Atom_type: integer_1 interger_2 ... to input script, for example Atom_type: 2 3" << endl;
            is_correct = false;
        }

        if(!is_correct)
        {
            cerr << "Flat_Membrane::validate_pore_input exiting" << endl;
            return -1;
        }
        return is_correct;
    }

    void is_Pore(Data& data, Atoms& mem)
    {
        cerr << endl;
        cerr << "Flat_Membrane::is_Pore" << endl;
        cerr << "- cell size = " << data.in.cell_size << endl;
        cerr << "- bead size = " << data.in.bead_size << endl;

        validate_pore_input(data);

        int number_of_cells = (data.in.sim_box.xhi - data.in.sim_box.xlo) / data.in.cell_size;
        double cell_size = (data.in.sim_box.xhi - data.in.sim_box.xlo) / number_of_cells;
        double inv_cell_size = 1.0 / cell_size;

        bool do_iterations = (data.in.cell_size < data.in.bead_size);
        int iterations = (0.5*data.in.bead_size / data.in.cell_size) +1;
        double R2 = 0.5*data.in.bead_size * 0.5*data.in.bead_size;

        bool is_pore = false;
        bool lattice[number_of_cells][number_of_cells];
        std::memset(lattice, false, sizeof(lattice)); // set all elements of latice to false -> false == no particle in cell

        Trajectory traj;
        Tensor_xyz box;
        traj.load(data.in.param["Trajectory_file"]);

        int cx=0, cy=0;
        int cell_ix=0, cell_iy=0;
        double cell_x_pos=0.0, cell_y_pos=0.0;
        double rx=0.0, ry=0.0;

        for(int i=0; i<traj.frame_count(); ++i)
        {
            mem.set_frame(traj[i]);
            box = traj.box_traj[i];
            cell_size = box.x / number_of_cells;
            inv_cell_size = number_of_cells / box.x; // box.x == box.y allways (enforced by lammps settings)

            std::memset(lattice, false, sizeof(lattice)); // set all elements of latice to false -> false == no particle in cell
            for(Atom a : mem)
            {
                for(int type : data.in.param_vector_int["Atom_type"])
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

            is_pore = false;
            for(int j=0; j<number_of_cells && !is_pore; ++j)
            {
                for(int k=0; k<number_of_cells; ++k)
                {
                    if( lattice[j][k] == false )
                    {
                        is_pore = true;
                        break;
                    }
                }
            }

            cout << i << " " << is_pore << endl; // true == 1

            /*if(true || i==1107) // debug
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
            }*/
        }
    }


    Lipids gen_flat_membrane(int num_lipids, int num_receptors, int mol_tag)
    {
        Lipids mem;

        int side_len = sqrt(num_lipids/2) +1;

        int count=0;
        double x=0.0,y=0.0;
        double z_up = 3.5, z_down=-3.5;

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

                    pos_up =   Tensor_xyz(x,y,z_up);
                    pos_down = Tensor_xyz(x,y,z_down);

                    mem.push_back(Lipid(pos_up,   dir_down, count,   mol_tag, Lipid::Leaflet::upper));
                    mem.push_back(Lipid(pos_down, dir_up,   count+1, mol_tag, Lipid::Leaflet::lower));
                    count+=2;
                }
            }
        }

        mem.convert_receptors(num_receptors);

        return mem;
    }

    void copy_Z(Data& data, Atoms& mem)
    {
        cerr << "Flat_Membrane::execute -> Copy_Z" << endl;

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
