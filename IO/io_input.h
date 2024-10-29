#ifndef IO_INPUT_H
#define IO_INPUT_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <cstdlib>
#include <algorithm>
#include <random>
#include <sstream>
#include <string.h>

#include <array>

#include "sim_box.h"
#include "atom.h"
#include "force_field.h"
#include "rng.h"

using namespace std;


enum class IO_Type { none, xyz, pdb, lammps_full };

std::ostream& operator<<(std::ostream& os, const IO_Type out)
{
    switch (out) {
        case IO_Type::xyz:
          os << "xyz"; return os;
        case IO_Type::pdb:
          os << "pdb"; return os;
        case IO_Type::lammps_full:
          os << "lammps_full"; return os;
        default:
          os << "none"; return os;
    }
    return os;
}




class IO{
public:
    IO() {}

    IO_Type type = IO_Type::none;

    void clear()
    {
        type = IO_Type::none;
    }

    friend std::istream& operator>>(std::istream& is, IO& out)
    {
        string str;
        is >> str;

        transform(str.begin(), str.end(), str.begin(), ::tolower);

        out.type = IO_Type::none;
        if (str == "xyz")
        {
            out.type = IO_Type::xyz;
        }
        if (str =="pdb")
        {
            out.type = IO_Type::pdb;
        }
        if (str == "lammps_full")
        {
            out.type = IO_Type::lammps_full;
        }

        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const IO& io)
    {
        os << io.type;
        return os;
    }
};



/**
 * @brief Populate
 * Defines the population of particles in the simulation box
 */
class Population
{
public:
	Population(){}

	bool random=false;
	int count=0;

	void clear()
	{
		random=false;
		count=0;
	}

    friend std::ostream& operator<<(std::ostream& os, const Population& pop)
    {
    	if(pop.random)
    	{
    		os << pop.count << " random";
    	}
    	else
    	{
    		os << "not defined";
    	}
        return os;
    }

	friend std::istream& operator>>(std::istream& is, Population& pop)
	{
	    string str;
	    is >> pop.count >> str;

	    transform(str.begin(), str.end(), str.begin(), ::tolower);

	    if (str == "random")
	    {
	    	pop.random = true;
	    }

	    return is;
	}
};




/**
 * @brief The Input class - Parameters for particle generation
 */
class IO_Input{
public:
    IO_Input() {}

    int i=0;

    // IO
    string infile; /// name of the filename with lammps_full atoms
    IO out; /// Output type - none, pdb, lammps_full, xyz
    IO in;  /// Input type  - none, pdb, lammps_full

    // particle identifier
    string gen_structure; /// keyword identifying the structure class

    // bead counts
    int num_of_beads=-1;
    int num_lig=-1;
    int subdiv_beads=-1;
    int subdiv_lig=-1;

    // atom types
    int chain_type=-1;
    int atom_type=1;

    // molecule types
    int mol_tag=-1;
    int mtag_1=-1;
    int mtag_2=-1;

    // patches
    vector<Atom> patches;
    Atom patch_1 = Atom(0,0,0,0);
    Atom patch_2 = Atom(0,0,0,0);

    // Ellipsoid, oblate spheroid
    myFloat b=0.0;
    myFloat c=0.0;

    // System data
    myFloat scale = 1.0;
    bool center=false;
    Tensor_xyz com_pos = Tensor_xyz(0.0, 0.0, 0.0);
    int seed=0;
    Tensor_xyz ivx = Tensor_xyz(0.0, 0.0, 0.0);
    bool fit = false;
    int offset = 1;

    // Population data
    Population population;

    // simulation box
    Simulation_Box sim_box;

    // Force-Field
    vector<LJ> bparam;
    vector<CosSQ> cparam;
    Force_Field ff;



    bool loadInput(string input, int i)
    {
        std::fstream fs( input, std::fstream::in );
        string line, what;
        stringstream ss;
        int len=0;

        this->i = i-1; // index of the file given to binary

        while( !fs.eof() ) // Lines in input
        {
        	what.clear();
            ss.flush();
            ss.clear();
            getline(fs, line);
            ss.str(line);

            ss >> what;

            // Load io
            if( what.compare("Load_file:") == 0 )         { ss >> infile; }
            if( what.compare("Output_type:") == 0 )       { ss >> out; }
            if( what.compare("Input_type:") == 0 )        { ss >> in; }

            // Load particle identifier
            if( what.compare("Particle_type:") == 0 )     { ss >> gen_structure; }

            // Load particle atom counts
            if( what.compare("Number_of_beads:") == 0 )   { ss >> num_of_beads; }
            if( what.compare("Number_of_ligands:") == 0 ) { ss >> num_lig; }
            if( what.compare("Subdiv_of_beads:") == 0 )   { ss >> subdiv_beads; }
            if( what.compare("Subdiv_of_ligands:") == 0 ) { ss >> subdiv_lig; }

            // Load atom and molecule types
            if( what.compare("Atom_type:") == 0 )         { ss >> atom_type; }
            if( what.compare("Mol_tag:") == 0 ) 		  { ss >> mol_tag; }
            if( what.compare("Chain_type:") == 0 ) 		  { ss >> chain_type; }

            // Load particle properties: aspect ration, patches
            if( what.compare("b:") == 0 )                 { ss >> b; }
            if( what.compare("c:") == 0 )                 { ss >> c; }
            if( what.compare("Patch:") == 0 )
            {
                patches.push_back(Atom());
                ss >> patches.back().pos.x >> patches.back().pos.y >> patches.back().pos.z;
                ss >> patches.back().vel.x >> patches.back().vel.y >> patches.back().vel.z;
                ss >> patches.back().type;
            }

            // Load system
            if( what.compare("Scale:") == 0 )             { ss >> scale; }
            if( what.compare("Center") == 0 )             { center=true; }
            if( what.compare("Position_shift:") == 0 )    { ss >> com_pos.x >> com_pos.y >> com_pos.z; }
            if( what.compare("Seed:") == 0 )              { ss >> seed; rng.seed(seed); }
            if( what.compare("Lammps_offset:") == 0 )     { ss >> offset; }

            // Load the simulation box
            if( what.compare("Sim_box:") == 0 )           { ss >> sim_box; }

            // Load the population of the particle
            if( what.compare("Populate:") == 0 )  		  { ss >> population; }

            // Load the force-field
            if( what.compare("ff_lj:") == 0 )             { LJ lj; ss >> lj; ff.lj[lj.type]=lj; }
            if( what.compare("ff_cos2:") == 0 )           { CosSQ cos; ss >> cos; ff.cos[cos.type]=cos; }

            // Load stuff for nanoparticle orientation and position
            if( what.compare("Align:") == 0 )  			{ ss >> mtag_1 >> mtag_2; }
            if( what.compare("Fit") == 0 )              { fit=true; }
            if( what.compare("Impact_vector:") == 0 )   { ss >> ivx.x >> ivx.y >> ivx.z; }
        }
        fs.close();

        return true;
    }

    string toString()
    {
        stringstream ss;

        if( !infile.empty() )
        {
            ss << "Load_file: " << infile << endl;
            ss << "Input_type:" << in << endl;
        }
        if( !gen_structure.empty() )
            ss << "Particle_type: " << gen_structure << endl;
        ss << "Output_type: " << out << endl;
        ss << "Number of beads: " << num_of_beads << endl;
        ss << "Number of ligands: " << num_lig << endl;
        ss << "Subdiv of beads: " << subdiv_beads << endl;
        ss << "Subdiv of ligands: " << subdiv_lig << endl;
        ss << "Scale: " << scale << endl;
        ss << "Offset: " << offset << endl;
        ss << "c: " << c << endl;
        ss << "Box: ( " << sim_box << " )" << endl;
        ss << "Position: (" << com_pos.x << ", " << com_pos.y << ", " << com_pos.z << ")" << endl;
        ss << "Patch_1: (" << patch_1.pos.x << "-" << patch_1.vel.x << ", " << patch_1.pos.y << "-" << patch_1.vel.y << ", " << patch_1.pos.z << "-" << patch_1.vel.z << ", " << patch_1.type << ")" << endl;
        ss << "Patch_1: (" << patch_2.pos.x << ", " << patch_2.pos.y << ", " << patch_2.pos.z << ", " << patch_2.type << ")" << endl;
        ss << "Populate: " << population << endl;

        ss << "Force-Field: " << ff << endl;

        return ss.str();
    }

    bool is_mtag_12()
    {
        if( mtag_1 == -1 || mtag_2 == -1 )
            return false;
        return true;
    }

    void clear()
    {
        // Persistent
        // out;
        // sim_box;
        // ff;

        gen_structure.clear();
        in.clear();

        num_of_beads=-1;
        num_lig=-1;
        subdiv_beads=-1;
        subdiv_lig=-1;

        scale=1.0;
        offset=1;
        b=0;
        c=0;

        atom_type=1;
        chain_type=-1;
        mol_tag=-1;

        center=false;
        fit = false;
        mtag_1=-1;
        mtag_2=-1;

        com_pos=Tensor_xyz(0.0, 0.0, 0.0);

        patches.clear();

        patch_1=Atom(0,0,0,0);
        patch_2=Atom(0,0,0,0);

        ivx=Tensor_xyz(0.0, 0.0, 0.0);

        population.clear();
        infile.clear();
        bparam.clear();
        cparam.clear();
    }

    void help()
    {
	    cout << "Lammps_offset: integer" << endl;
	    cout << " - offset the generated structure for manual insertion into another lammps structure file" << endl;
	    cout << "Load_file: filename" << endl;


	    cout << "Num_of_beads: integer" << endl;
	    cout << " - Particle_type: 1,6,7 = number of beads per edge" << endl;
	    cout << " - other Particle_type = number of beads for entire nanoparticle/structure" << endl;
	    cout << "Scale: floating_point_number - nanoparticle radius" << endl;
	    cout << "c: float " << endl;
	    cout << " - Particle_type: 4 = oblate spheroid < 1, prolate spheroid > 1, 1.0 - ERROR, not defined" << endl;
	    cout << " - Particle_type: 3 = width of patch (0.0 to 2.0)" << endl;
	    cout << "Number_of_ligands: integer" << endl;
	    cout << "Mol_tag: integer" << endl;
	    cout << " - change mol_tag of generated/loaded structure " << endl;
	    cout << "Atom_type: integer" << endl;
	    cout << " - Atom_type of generated structure, if structure has more atom_types they are incremented from provided value " << endl;
	    cout << "Janus: float float float" << endl;

	    cout << "\nPosition/Box properties category:" << endl;
	    cout << "Box: float float float float float float" << endl;
	    cout << "position_shift: float float float" << endl;
	    cout << "Center = centers particles at 0.0" << endl;
	    cout << "Align: mol_tag_1 mol_tag_2 integer" << endl;
	    cout << " - align mol_tag_1 with x-axis" << endl;
	    cout << " - center mol_tag_2 around z axis" << endl;
	    cout << "Impact_vector: float float float" << endl;
	    cout << "Fit" << endl;
	    cout << " - positions the loaded/generated structure next to previosly generated/loadedd structure" << endl;
	    cout << " - used for ideal collision position of two liposomes and a nanoparticle" << endl;

	    cout << "\nForce-Field category:" << endl;
	    cout << "Beads_lj/cut:" << endl;

	    cout << "\nSeed: integer = random generator" << endl;
    }
};


#endif // IO_INPUT_H
