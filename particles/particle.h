#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "xdrfile-1.1.4/include/xdrfile.h"
#include "xdrfile-1.1.4/include/xdrfile_xtc.h"

#include<stdio.h>

#include "data.h"
#include "atom.h"


using namespace std;




/**
 * @brief The Particle class - Abstract class of Particle - a assembly made from several beads
 */
class Particle{
public:
    inline static const string keyword = "Abstract class of Particle";
    const string name = "Abstract class of Particle";

    Simulation_Box box;

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
        scale( data.in.scale );
        move( data.in.com_pos );
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

        data.coll_beads.push_back(beads);
        data.coll_bonds.push_back(bonds);
        data.coll_angles.push_back(angles);

        data.all_sigma = this->sigma;
        data.all_sigma_size = this->sigma_size;
        data.all_sigma_cosatt = this->sigma_cosatt;

        beads.clear();
        bonds.clear();
        angles.clear();
    }

    virtual string help()
    {
        stringstream ss;
        ss << "Abstract class Particle\n";
        ss << "Contains functions intended for inheritance\n";
        ss << "Does not generate anything\n";
        return ss.str();
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

    			temp.set_mol_tag(copy[0].mol_tag+i);

    			// rotate copied structure
    			Atoms::clusterRotate_random(temp, 180.0*degToRad);

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
    	cerr << "move " << move.pos.x << " " << move.pos.y << " " << move.pos.z << endl;
    }

    void scale(double scale)
    {
    	beads.scale(scale);
    	cerr << "scale " << scale << endl;
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
        beads.push_back(Atom(0,0,0,data.in.atom_type[0],data.in.mol_tag));
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
