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
class Particles{
public:
    inline static const string keyword = "Abstract class of Particle";
    const string name = "Abstract class of Particle";

    Simluation_Box box;

    const double degToRad = 0.0174532925;

    /// Local data for structure, need to save to Data class for output
    Atoms beads;
    Bonds bonds;
    vector<Angle> angles;
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
    Particles() {}
    Particles(string name) : name(name) {}
    virtual ~Particles() {}


    /**
     * @brief generate - Generate particle
     * - Particle need to be generated in scale and in position (data.in.scale and data.in.com_pos)
     * - Position overlap check to already existing particles in data.all_beads
     *
     * @param data
     */
    virtual void generate( Data& data )=0;

    void populate(Data& data)
    {
    	Atoms copy(beads);
    	beads.clear();
    	Atoms temp;
    	Atom com_pos;

    	if( data.in.population.random )
    	{
    		for(int i=0; i < data.in.population.count; ++i)
    		{
    			com_pos = data.in.sim_box.get_random_pos();
    			temp = copy;
    			temp.move(com_pos);

    			int tries=0;
    			while( beads.is_overlap(temp) )
    			{
    				com_pos = data.in.sim_box.get_random_pos();
    				temp = copy;
    				temp.move(com_pos);

    				if(tries > 1000)
    				{
    					cerr << "Can't generate any more particles, simulation box is full, tries " << tries << endl;
    					exit(3);
    				}
    				++tries;
    			}

    			temp.set_mol_tag(i+1);
    			beads.insert( beads.end(), temp.begin(), temp.end() );
    			cerr << "Generating population, particle " << i << " of " << data.in.population.count << ", successful on try " << tries << endl;
    		}
    	}
    	cerr << "populate " << data.in.population << endl;
    }

    void move(Atom move)
    {
    	beads.move(move);
    	cerr << "move " << move.x << " " << move.y << " " << move.z << endl;
    }

    void scale(double scale)
    {
    	beads.scale(scale);
    	cerr << "scale " << scale << endl;
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

    void add(Data& data)
    {
        bool overlap = false;
        for(auto& item : this->beads)
        {
            if( data.isOverlap(item) )
            {
                overlap = true;
            }
        }

        data.all_beads.insert(data.all_beads.end(), this->beads.begin(), this->beads.end());
        data.all_bonds.insert(data.all_bonds.end(), this->bonds.begin(), this->bonds.end());
        data.all_angles.insert(data.all_angles.end(), this->angles.begin(), this->angles.end());

        data.all_sigma = this->sigma;
        data.all_sigma_size = this->sigma_size;
        data.all_sigma_cosatt = this->sigma_cosatt;

        if(overlap)
        {
            cerr << "WARNING: overlap detected" << endl;
        }
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
};


class Empty_Particle : public Particles
{
public:
    inline static const string keyword = "empty";
    const string name = "empty";
    Empty_Particle() : Particles("empty") {}
    void generate( Data& data ) {}

    string help()
    {
    	stringstream ss;

    	ss << "Particle_type: empty\n";
    	ss << "Output_type: pdb # other keywords: xyz pdb lammps_full\n";

    	return ss.str();
    }
};

#endif // PARTICLE_H
