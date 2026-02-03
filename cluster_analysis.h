#ifndef CLUSTER_ANALYSIS_H
#define CLUSTER_ANALYSIS_H

#include <vector>
#include "atom.h"
#include "lipid.h"




class Unassigned : public vector<int>
{
public:
    Unassigned(Atoms& a)
    {
        unassigned_size = a.size();
        resize(a.size());
        for(int i=0; i<size(); ++i)
        {
            at(i) = i;
        }
    }

    Unassigned(Atoms& a, bool only_tails)
    {
        only_tails = true;
        Lipid test;
        resize(a.size());

        for(int i=0; i<a.size(); ++i)
        {
            if( test.is_tail(a[i].type) )
            {
                at(unassigned_size) = i;
                ++unassigned_size;
            }
        }
    }

    bool only_tails = false;
    int unassigned_size = 0;

    void remove(int index)
    {
        at(index) = at(unassigned_size-1); // move last id to front
        --unassigned_size; // decrease size by 1
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
    Clusters(){}

    int analyze(Atoms& a, double cutoff)
    {
        bool only_tails = true;
        Unassigned un_id(a, only_tails); // list of particles not assigned to a cluster
        int safety = 10*1000;
        int count=0;

        while(un_id.unassigned_size > 0 && count<safety )
        {
            push_back(Cluster( un_id.at(0) )); // add the first particle
            un_id.remove(0);
            back().fill_cluster(a, un_id, cutoff);
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
