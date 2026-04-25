#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <sstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <string>

#include "data.h"
#include "atom.h"
#include "cell_list.h"


using namespace std;




/**
 * @brief The Particle class - Abstract class of Particle - a assembly made from several beads
 */
class Particle{
public:
    inline static const string keyword = "Abstract class of Particle";
    const string name = "Abstract class of Particle";

    Simulation_Box box;

    Cell_List cell_list;
    Cell_List coll_cell_list;

    const double degToRad = 0.0174532925;

    /// Local data for structure, need to save to Data class for output
    Atoms beads;
    Bonds bonds;
    Angles angles;
    vector<LJ> bparam;

    int typeNano = 1;
    int typeLig = 2;
    const int typeTemp = 99;

    const int orientX = 101;
    const int orientY = 102;
    const int orientZ = 103;

    /// Force-field Data
    int sigma_size = 0;
    array<array<double, 100>, 100> sigma;   
    array< array<bool, 100>, 100> sigma_cosatt;

    //
    // Class stuff
    //
    Particle() {}
    Particle(string name) : name(name) {}
    virtual ~Particle() {}


    /**
     * @brief generate - Generate particle
     * - Particle need to be generated in scale and in position (data.in.scale and data.in.com_pos)
     * - Position overlap check to already existing particles in data.all_beads
     *
     * @param data
     */
    virtual void generate( Data& data )=0;

    void modify( Data& data )
    {
        cerr << "Particle::modify" << endl;
        if( data.in.p_float.contains("Scale") )          { scale( data.in.p_float["Scale"] ); }
        if( data.in.p_tensor.contains("Position_shift")) { move(  data.in.p_tensor["Position_shift"] ); }
    }

    void make_persistent(Data& data)
    {
        for(Atoms& bds : data.coll_beads)
        {
            if(bds.is_overlap(beads, data.in.ff))
            {
                cerr << "WARNING: overlap detected" << endl;
            }
        }

        add_beads_to_coll(data);
        data.coll_bonds.push_back(bonds);
        data.coll_angles.push_back(angles);

        data.all_sigma = this->sigma;
        data.all_sigma_size = this->sigma_size;
        data.all_sigma_cosatt = this->sigma_cosatt;

        beads.clear();
        bonds.clear();
        angles.clear();
    }

    void add_beads_to_coll(Data& data)
    {
        beads.offset_N(get_coll_N(data));

        if(data.in.p_int.contains("ID") && data.id_map.contains( data.in.p_int["ID"] ) && !data.coll_beads.empty()) // adding to existing collection
        {
            // Happens when: 1. defined ID, 2. loaded (even empty) file, 3. generated a particle
            if(data.id_map.contains( data.in.p_int["ID"] ))
            {
                int sys_id = data.id_map[ data.in.p_int["ID"] ];
                if(data.coll_beads.size() > sys_id)
                {
                    data.coll_beads[sys_id].insert(data.coll_beads[sys_id].end(), beads.begin(), beads.end());
                }
                else
                {
                    cerr << "ERROR, Particle::add_beads_to_coll data.coll_beads does not contain sys_id index" << endl;
                    exit(1);
                }
            }
            else
            {
                cerr << "ID defined, but not stored in id_map" << endl;
                exit(1);
            }
        }
        else // adding a new collection
        {
            data.coll_beads.push_back(beads);

            if(data.in.p_int.contains("ID")) // ID defined, but not yet used (by data loading)
            {
                data.id_map[ data.in.p_int["ID"] ] = data.coll_beads.size()-1;
            }
        }
    }

    int get_coll_N(Data& data)
    {
        int offset_N = 0;

        if(!data.coll_beads.empty())
        {
            for(Atoms& coll : data.coll_beads)
            {
                if(!coll.empty() && coll.back().N > offset_N)
                {
                    offset_N = coll.back().N;
                }
            }
        }

        return offset_N;
    }

    virtual string help()
    {
        return "";
    }

    void printSigma()
    {
        cerr << "printSigma:" << endl;
        for(unsigned int i=0; i<sigma_size; ++i)
        {
            cerr << "[" << i+1 << "][" << i+1 << "] = "  << sigma[i][i] << endl;
        }
    }

    void printSigma(int j)
    {
        cerr << "printSigma type " << j << ":" << endl;
        for(unsigned int i=0; i<sigma_size; ++i)
        {
            cerr << "[" << j+1 << "][" << i+1 << "] = "  << sigma[j][i] << " == " << sigma[i][j] << endl;
        }
    }

    bool isSame (vector<Atom>& container, Atom& push)
    {
        bool same = false;

        for(unsigned int q=0; q< container.size(); q++) {
            if(container[q].isAproxSame(push)) {
                same = true;
            }
        }

        return !same;
    }

    virtual void populate(Data& data)
    {
    	if( data.in.population.random )
    	{
            //
            // Report
            //
            if(beads.empty())
            {
                cerr << "No particle created, can't populate system" << endl;
                exit(-1);
            }
            else
            {
                cerr << "Copying particle of " << beads.size() << " atoms " << data.in.population.count << " times" << endl;
            }

            Atoms copy(beads);
            beads.clear();
            Atoms temp;
            Tensor_xyz com_pos;

            //
            // Copy, move, rotate
            //
    		for(int i=0; i < data.in.population.count; ++i)
    		{
    			com_pos = data.in.sim_box.get_random_pos();
    			temp = copy;
    			temp.move(com_pos);

    			int tries=0;
                while( beads.is_overlap(temp, data.in.ff) || data.is_overlap(temp, data.in.ff) )
    			{
    				com_pos = data.in.sim_box.get_random_pos();
    				temp = copy;

                    temp.move(com_pos);

                    if(tries > 1000)
    				{
                        cerr << "Particles::populate -> Can't generate any more particles, simulation box is full, tries" << tries << endl;
                        exit(1);
    				}
    				++tries;
    			}

                // increment N
                if(beads.empty())
                {
                    int count = 1; // N=1 by default
                    for(Atom& a : temp)
                    {
                        a.N = count;
                        ++count;
                    }
                }
                else
                {
                    for(Atom& a : temp)
                    {
                        a.N += beads.back().N;
                    }
                }

                // mol tag increment
                //temp.set_mol_tag(copy[0].mol_tag+i);

    			// rotate copied structure
                //Atoms::clusterRotate_random(temp, 180.0*degToRad);

    			beads.insert( beads.end(), temp.begin(), temp.end() );

    			cerr << "Generating population, particle " << i << " of " << data.in.population.count << ", successful on try " << tries << endl;
    		}
    	}
    	cerr << "populate " << data.in.population << endl;
    }

private:

    void move(Atom move)
    {
    	beads.move(move);
        cerr << "move: " << move.pos << endl;
    }

    void scale(double scale)
    {
    	beads.scale(scale);
        cerr << "scale: " << scale << endl;
    }
};


class Empty_Particle : public Particle
{
public:
    inline static const string keyword = "empty";
    const string name = "empty";
    Empty_Particle() : Particle("empty") {}
    void generate( Data& data ) {}

    string help()
    {
    	stringstream ss;

    	ss << "Particle_type: empty\n";
    	ss << "Output_type: pdb # other keywords: xyz pdb lammps_full\n";

    	return ss.str();
    }
};

class Monomer : public Particle
{
public:
    inline static const string keyword = "monomer";
    const string name = "monomer";
    Monomer() : Particle("monomer") {}

    void generate( Data& data )
    {
        beads.push_back( Atom(0.0, 0.0, 0.0, data.in.p_vec_int["Atom_type"][0], data.in.p_int["Mol_tag"]) );
    }

    string help()
    {
    	stringstream ss;

    	ss << "Particle_type: monomer\n";
    	ss << "Output_type: pdb # other keywords: xyz pdb lammps_full\n";

    	return ss.str();
    }
};

#endif // PARTICLE_H
