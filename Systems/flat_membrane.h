#ifndef FLAT_MEMBRANE_H
#define FLAT_MEMBRANE_H

#include "system_base.h"
#include "atom.h"

#include "lipid.h"

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
            //
            // Analysis only on lipid tails
            //
            Lipid test;

            cell_list.init(data);
            cell_list.add(mem);



            vector<int> unassigned; // list of particles not assigned to a cluster
            int unassigned_size = mem.size();
            unassigned.resize(mem.size());
            for(int i=0; i<unassigned.size(); ++i)
            {
                unassigned[i] = i;
            }

            vector<int> cluster; // list of particles of cluster
            cluster.push_back(0); // add the first particle
            double cutoff=1.6;
            int particle_ID = 0;
            bool added=false;
            for(int i=0; i<cluster.size(); ++i) // cluster.size() dynamically evaluated each time
            {
                particle_ID = cluster[i];
                for(int j=1; j<mem.size(); ++j)
                {
                    if(mem[particle_ID].distSQ(mem[j]) < cutoff*cutoff)
                    {
                        added = false;
                        for(int k=0; k<cluster.size(); ++k)
                        {
                            if(cluster[k] == j)
                            {
                                added = true;
                            }
                        }
                        if(!added)
                        {
                            cluster.push_back(j);
                            //cerr << "Push " << particle_ID << "::" << j << " " << mem[particle_ID].dist(mem[j]) << endl;
                        }
                    }
                }
            }
            cerr << "Cluster size: " << cluster.size() << endl;
            exit(1);
        }
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
