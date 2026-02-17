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

    int unassigned_size = 0;

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

    void fill_cluster(Atoms& a, Unassigned& id, double cutoff)
    {
        int particle_ID = 0;
        bool added=false;
        for(int i=0; i<size(); ++i) // cluster.size() dynamically evaluated each time
        {
            particle_ID = at(i);
            for(int j=0; j<id.unassigned_size; ++j) // unassigned particles
            {
                if(a[particle_ID].distSQ( a[id.at(j)] ) < cutoff*cutoff)
                {
                    added = false;
                    for(int k=0; k<size(); ++k)
                    {
                        if(at(k) == id.at(j))
                        {
                            added = true;
                        }
                    }
                    if(!added)
                    {
                        push_back(id.at(j));
                        id.remove(j);
                        //cerr << "Push " << particle_ID << "::" << j << " " << mem[particle_ID].dist(mem[j]) << endl;
                    }
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

    Unassigned un_id; // list of particles not assigned to a cluster of a given type

    int analyze(Atoms& a, double cutoff)
    {
        int safety = 10*1000;
        int count=0;

        Unassigned un_id_temp = un_id;

        while(un_id_temp.unassigned_size > 0 && count<safety )
        {
            push_back(Cluster( un_id_temp.at(0) )); // add the first particle
            un_id_temp.remove(0);
            back().fill_cluster(a, un_id_temp, cutoff);
            count++;
        }

        int total_particle_count = 0;
        int cluster_count = size();

        for(int i=0; i<size(); ++i)
        {
            total_particle_count += at(i).size();
        }

        /*if(a.size() != total_particle_count) // does not take into account tail selection
        {
            cerr << "Total particle count = " << a.size() << " != particle count in clusters " << total_particle_count << endl;
        }*/

        return cluster_count;
    }
};

#endif // CLUSTER_ANALYSIS_H
