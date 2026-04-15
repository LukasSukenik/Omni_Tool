#ifndef CLUSTER_ANALYSIS_H
#define CLUSTER_ANALYSIS_H

#include <vector>
#include "atom.h"
#include "lipid.h"
#include "cell_list.h"
#include "sim_box.h"

class Cluster_Basic; // declaration of classes, must be done for using Clusters
class Cluster_Cell_List;

class Clusters_Basic;
class Clusters_Cell_List;

using Cluster = Cluster_Cell_List;
using Clusters = Clusters_Cell_List;


//
// Interface - TODO
//
class Clusters_interface
{
public:
    virtual int analyze(Atoms& topo, Simulation_Box& sim_box, double cutoff) =0;
    virtual size_t cluster_count() =0;
};

//
// Cluster analysis - accelerated with cell list
//
class In_Cluster : public vector<bool>
{
public:
    In_Cluster(Atoms& topo, vector<int> types)
    {
        resize(topo.size(), true);
        for(size_t i=0; i<topo.size(); ++i)
        {
            for(size_t j=0; j<types.size(); ++j)
            {
                if(topo[i].type == types[j])
                {
                    this->at(i) = false;
                }
            }
        }
    }
};

class Cluster_Cell_List : public vector<int>
{
public:
    bool is_filled = false;

    Cluster_Cell_List(){}

    void add_particle(int atom_ID) { push_back(atom_ID); }

    void fill_cluster(Atoms& topo, Cell_List& cell_list, In_Cluster& in_cluster, double cutoff)
    {
        int c_p_ID = 0; // cluster particle ID

        for(int i=0; i<size(); ++i) // Loop over cluster particles, size() dynamically evaluated
        {
            c_p_ID = at(i);
            cell_list.set_neighbors_pbc( topo[c_p_ID].pos );
            for(int j=0; j<cell_list.neighbors.size(); ++j) // Loop over cells
            {
                for(int nei_ID : (*cell_list.neighbors[j]) ) // Loop over neigbor particles in those cells
                {
                    if(in_cluster[ nei_ID ] == false && topo[ c_p_ID ].distSQ_pbc( topo[ nei_ID ], cell_list.pbc, cell_list.pbc_inv ) < cutoff*cutoff)
                    {
                        push_back( nei_ID );
                        in_cluster[ nei_ID ] = true;
                    }
                }
            }
        }
        is_filled = true;
    }

private:
    void fill_cluster_without_cell_list(Atoms& topo, In_Cluster& in_cluster, double cutoff)
    {
        int cluster_part_ID = 0;
        for(int i=0; i<size(); ++i) // Loop over cluster particles, size() dynamically evaluated
        {
            cluster_part_ID = at(i);
            for(int j=0; j<topo.size(); ++j) // Loop over remaining un-assigned particles
            {
                if(in_cluster[j] == false && topo[cluster_part_ID].distSQ( topo[ j ] ) < cutoff*cutoff)
                {
                    push_back(j);
                    in_cluster[j] = true;
                }
            }
        }
        is_filled = true;
    }
};

class Clusters_Cell_List : public vector< Cluster_Cell_List >
{
public:
    vector<int> types;

    Clusters_Cell_List(){}
    Clusters_Cell_List(Atoms& topo, vector<int> types) : types(types) {}

    int analyze(Atoms& topo, Simulation_Box& sim_box, double cutoff)
    {
        In_Cluster in_cluster(topo, types); // set all to true, then for types vector set to false
        Cell_List cell_list;
        cell_list.init( sim_box ); // allocate memory
        cell_list.add( topo ); // populate the cell list

        for(size_t i=0; i<topo.size(); ++i)
        {
            if(in_cluster[i] == false)
            {
                push_back(Cluster_Cell_List());
                back().add_particle(i);
                in_cluster[i] = true;
                back().fill_cluster(topo, cell_list, in_cluster, cutoff);
            }
        }

        cell_list.delete_all();

        return cluster_count();
    }

    size_t cluster_count() { return size(); }
};





//
// Deprecated cluster analysis classes - slow performance
//

/**
 * @brief The Unassigned class - Helper class to class Clusters
 * - keep track of atoms unassigned to a clusters
 */
class Unassigned_int : public vector<int>
{
private:
    int unassigned_size = 0;
public:
    Unassigned_int(Atoms& topo) // O(n) complexity, all types
    {
        unassigned_size = topo.size();
        resize(topo.size());
        for(int i=0; i<size(); ++i)
        {
            at(i) = i;
        }
    }

    Unassigned_int(Atoms& a, vector<int> types) // O(n) complexity
    {
        Lipid test;
        resize(a.size());

        for(int i=0; i<a.size(); ++i)
        {
            for(int type : types)
            {
                if( type == a[i].type )
                {
                    at(unassigned_size) = i;
                    ++unassigned_size;
                }
            }
        }
    }

    int unassigned_count()
    {
        return unassigned_size;
    }

    void remove(int index)
    {
        at(index) = at(unassigned_size-1); // move last id to front
        --unassigned_size; // decrease size by 1
    }

    // Copy constructor
    Unassigned_int(const Unassigned_int& o) : std::vector<int>(o), unassigned_size(o.unassigned_size) {}

    // Copy assignment operator
    Unassigned_int& operator=(const Unassigned_int& o)
    {
        if (this != &o)
        {
            std::vector<int>::operator=(o);
            this->unassigned_size = o.unassigned_size;
        }
        return *this;
    }
};




class Cluster_Basic : public vector<int>
{
public:
    Cluster_Basic(){}
    Cluster_Basic(int atom_ID)
    {
        push_back(atom_ID);
    }

    bool is_filled = false;

    void fill_cluster(Atoms& topo, Unassigned_int& unassigner_id, double cutoff)
    {
        int cluster_part_ID = 0;
        for(int i=0; i<size(); ++i) // Loop over cluster particles, cluster.size() dynamically evaluated
        {
            cluster_part_ID = at(i);
            for(int j=0; j<unassigner_id.unassigned_count(); ++j) // Loop over remaining un-assigned particles
            {
                if(topo[cluster_part_ID].distSQ( topo[ unassigner_id[j] ] ) < cutoff*cutoff)
                {
                    push_back(unassigner_id.at(j));
                    unassigner_id.remove(j);
                }
            }
        }
        is_filled = true;
    }
};




class Clusters_Basic : public vector< Cluster_Basic >
{
public:
    Clusters_Basic(Atoms& a, vector<int> types) : un_id(Unassigned_int(a, types)) {}

    const Unassigned_int un_id; // list of all particles, algorithm requires them not assigned to a cluster

    size_t cluster_count()
    {
        return size();
    }

    /**
     * @brief analyze - returns the number of clusters in topo
     * @param topo
     * @param cutoff
     * @return
     */
    int analyze(Atoms& topo, Simulation_Box& sim_box, double cutoff)
    {
        size_t max_cluster_count = 10*1000;
        Unassigned_int un_id_temp = un_id;

        while(un_id_temp.unassigned_count() > 0 && cluster_count()<max_cluster_count ) // loop over clusters
        {
            push_back(Cluster_Basic( un_id_temp.at(0) )); // new cluster with 1 unassigned particle
            un_id_temp.remove(0); // remove unassined particle
            back().fill_cluster(topo, un_id_temp, cutoff); // back() is the new cluster, fill is O(N^3)
        }
        return cluster_count();
    }
};

#endif // CLUSTER_ANALYSIS_H
