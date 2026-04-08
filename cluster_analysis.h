#ifndef CLUSTER_ANALYSIS_H
#define CLUSTER_ANALYSIS_H

#include <vector>
#include "atom.h"
#include "lipid.h"



/**
 * @brief The Unassigned class - Helper class to class Clusters
 * - keep track of atoms unassigned to a clusters
 */
class Unassigned : public vector<int>
{
private:
    int unassigned_size = 0;
public:
    Unassigned(Atoms& a) // O(n) complexity, all types
    {
        unassigned_size = a.size();
        resize(a.size());
        for(int i=0; i<size(); ++i)
        {
            at(i) = i;
        }
    }

    Unassigned(Atoms& a, vector<int> types) // O(n) complexity
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
    Unassigned(const Unassigned& o) : std::vector<int>(o), unassigned_size(o.unassigned_size) {}

    // Copy assignment operator
    Unassigned& operator=(const Unassigned& o)
    {
        if (this != &o)
        {
            std::vector<int>::operator=(o);
            this->unassigned_size = o.unassigned_size;
        }
        return *this;
    }
};




class Cluster : public vector<int>
{
public:
    Cluster(){}
    Cluster(int atom_ID)
    {
        push_back(atom_ID);
    }

    bool is_filled = false;

    void fill_cluster(Atoms& topo, Unassigned& unassigner_id, double cutoff)
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




class Clusters : public vector< Cluster >
{
public:
    Clusters(Atoms& a, vector<int> types) : un_id(Unassigned(a, types)) {}

    const Unassigned un_id; // list of all particles, algorithm requires them not assigned to a cluster

    size_t cluster_count()
    {
        return size();
    }

    int analyze(Atoms& topo, double cutoff)
    {
        size_t max_cluster_count = 10*1000;
        Unassigned un_id_temp = un_id;

        while(un_id_temp.unassigned_count() > 0 && cluster_count()<max_cluster_count ) // loop over clusters
        {
            push_back(Cluster( un_id_temp.at(0) )); // new cluster with 1 unassigned particle
            un_id_temp.remove(0); // remove unassined particle
            back().fill_cluster(topo, un_id_temp, cutoff); // back() is the new cluster, fill is O(N^3)
        }
        return cluster_count();
    }
};

#endif // CLUSTER_ANALYSIS_H
