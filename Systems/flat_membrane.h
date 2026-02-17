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

    void execute(Data& data)
    {
        int sys_id = data.id_map[data.in.id];
        Atoms& mem = data.coll_beads[sys_id];
        vector<int> mol_tags = mem.get_Mol_Types();

        if(data.in.system_function.compare("Copy_Z") == 0) { copy_Z(data, mem); }
        if(data.in.system_function.compare("Cluster_Analysis") == 0)
        {
            cerr << endl;
            cerr << "Cluster analysis: types ";

            for(int type : data.in.cluster_types)
            {
                cerr << type << ", ";
            }
            cerr << " | cutoff: " << data.in.cluster_cutoff << endl;;

            Clusters clusters(mem, data.in.cluster_types); // list of particle indexes
            Trajectory traj;

            traj.load(data.in.trajectory);

            for(int i=0; i<traj.frame_count(); ++i)
            {
                mem.set_frame(traj[i]);
                cout << i << " " << clusters.analyze(mem, data.in.cluster_cutoff) << endl;
                clusters.clear();
            }

        }
        if(data.in.system_function.compare("Calc_Z_Dist") == 0) { calc_Z_Dist(data, mem); }
        if(data.in.system_function.compare("Calc_Pore") == 0) { is_Pore(data, mem); }
    }


    void generate( Data& data )
    {
        Lipids membrane = gen_flat_membrane(data.in.num_lipids, data.in.num_rec, data.in.mol_tag);

        for(Lipid& lip : membrane)
        {
            beads.insert(beads.end(), lip.part.begin(), lip.part.end());
            bonds.insert(bonds.end(), lip.bond.begin(), lip.bond.end());
        }

        cerr << "suggesting box size: " << - 0.5*(sqrt(data.in.num_lipids/2)+1) << " to " << 0.5*(sqrt(data.in.num_lipids/2)+1) << endl;
    }

private:
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
        traj.load(data.in.trajectory);

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




    void is_Pore(Data& data, Atoms& mem)
    {
        //
        // Calculate whether the stalk has pore
        //
        // Project lipid tails along Z axis to a boolean 2D array
        // - when each element is set to true == there is a bead projected to the element
        // - when the element is set to false == pore
        //
        // The boolean 2d Array is mapped to the simulation box in the z_normal plane
        // - the size of each element - a square is set
        //

        // Accessing box size in trajectory
        Trajectory traj;
        traj.load(data.in.trajectory);
        cerr << traj.box_traj[0].x << endl;

        // This is similar to the histogram function in logic
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

        for(int i=0; i<mem.size()/4; ++i)
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
        for(int i=0; i<a.size(); ++i)
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
