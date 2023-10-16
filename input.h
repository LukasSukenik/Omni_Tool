#ifndef INPUT_H
#define INPUT_H

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

#include "atom.h"
#include "force_field.h"

using namespace std;




class Simluation_Box
{
public:
	Simluation_Box() {}

    myFloat xlo = 0.0;
    myFloat xhi = 0.0;
    myFloat ylo = 0.0;
    myFloat yhi = 0.0;
    myFloat zlo = 0.0;
    myFloat zhi = 0.0;

    inline double volume() {
        return (xhi-xlo)*(yhi-ylo)*(zhi-zlo);
    }

    Atom get_random_pos()
    {
    	Atom a;
    	a.x = ran() * (xhi-xlo) + xlo;
    	a.y = ran() * (yhi-ylo) + ylo;
    	a.z = ran() * (zhi-zlo) + zlo;
    	return a;
    }

    Atom usePBC(Atom& pos_orig, double scale) const {
        Atom pos = pos_orig;

        while (pos.x < 0.0) {
            pos.x += xhi/scale;
        }
        while (pos.x > xhi/scale) {
            pos.x -= xhi/scale;
        }
        while (pos.y < 0.0) {
            pos.y += yhi/scale;
        }
        while (pos.y > yhi/scale) {
            pos.y -= yhi/scale;
        }
        while (pos.z < 0.0) {
            pos.z += zhi/scale;
        }
        while (pos.z > zhi/scale) {
            pos.z -= zhi/scale;
        }
        return pos;
    }

    friend std::istream& operator>>(std::istream& is, Simluation_Box& sim_box)
    {
    	is >> sim_box.xlo >> sim_box.xhi >> sim_box.ylo >> sim_box.yhi >> sim_box.zlo >> sim_box.zhi;
        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const Simluation_Box& box)
    {
    	os << box.xlo << ", " << box.xhi << ", " << box.ylo << ", " << box.yhi << ", " << box.zlo << ", " << box.zhi;
        return os;
    }
};




enum class Output_Type { none, xyz, pdb, lammps_full };

std::ostream& operator<<(std::ostream& os, const Output_Type out)
{
    switch (out) {
        case Output_Type::xyz:
          os << "xyz"; return os;
        case Output_Type::pdb:
          os << "pdb"; return os;
        case Output_Type::lammps_full:
          os << "lammps_full"; return os;
        default:
          os << "none"; return os;
    }
    return os;
}




class Output{
public:
    Output() {}

    Output_Type type;

    void clear()
    {
        type = Output_Type::none;
    }

    friend std::istream& operator>>(std::istream& is, Output& out)
    {
        string str;
        is >> str;

        transform(str.begin(), str.end(), str.begin(), ::tolower);

        out.type = Output_Type::none;
        if (str == "xyz")
        {
            out.type = Output_Type::xyz;
        }
        if (str =="pdb")
        {
            out.type = Output_Type::pdb;
        }
        if (str == "lammps_full")
        {
            out.type = Output_Type::lammps_full;
        }

        return is;
    }

    friend std::ostream& operator<<(std::ostream& os, const Output& out)
    {
        os << out.type;
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
class Input{
public:
    Input() {}

    //
    // Data updated per input file
    //
    Output out; /// Output type - none, pdb, lammps_full, xyz
    string gen_structure; /// keyword identifying the structure class
    string infile; /// name of the filename with lammps_full atoms

    /// Transformation of the generated structure
    myFloat scale = 0.0;
    bool center=false;
    Atom com_pos = Atom(0.0, 0.0, 0.0);

    // structure variables
    int num_of_beads=-1;
    int num_lig=-1;
    myFloat c;
    vector<Atom> patches;
    Atom patch_1 = Atom(0,0,0,0);
    Atom patch_2 = Atom(0,0,0,0);

    /// Atom Types
    int chain_type=-1;
    int mol_tag=-1;
    int atom_type=1;
    int mtag_1=-1;
    int mtag_2=-1;

    /// Population data
    Population population;

    int seed=0;

    //
    // Persistent data
    //
    Simluation_Box sim_box;

    Atom ivx = Atom(0.0, 0.0, 0.0);
    bool fit = false;

    int offset = 1;

    // Force-Field
    Force_Field ff;
    vector<LJ> bparam;
    vector<CosSQ> cparam;

    bool loadInput(string input)
    {
        std::fstream fs( input, std::fstream::in );
        string line, what;
        stringstream ss;
        int len=0;

        while( !fs.eof() ) // Lines in input
        {
            ss.flush();
            ss.clear();
            getline(fs, line);
            ss.str(line);

            ss >> what;
            if( what.compare("Output_type:") == 0 )     { ss >> out; }
            if( what.compare("Particle_type:") == 0 )   { ss >> gen_structure; }
            if( what.compare("Load_file:") == 0 )  { ss >> infile; }

            if( what.compare("Scale:") == 0 )           { ss >> scale; }
            if( what.compare("Center") == 0 )  { center=true; }
            if( what.compare("Position_shift:") == 0 )  { ss >> com_pos.x >> com_pos.y >> com_pos.z; }

            if( what.compare("Num_of_beads:") == 0 )    { ss >> num_of_beads; }
            if( what.compare("Number_of_ligands:") == 0 ) { ss >> num_lig; }
            if( what.compare("c:") == 0 )               { ss >> c; }
            if( what.compare("Patch:") == 0 )  		{ patches.push_back(Atom()); ss >> patches.back().x >> patches.back().y >> patches.back().z >> patches.back().vx >> patches.back().vy >> patches.back().vz >> patches.back().type; }


            if( what.compare("Chain_type:") == 0 ) 		{ ss >> chain_type; }
            if( what.compare("Mol_tag:") == 0 ) 		{ ss >> mol_tag; }
            if( what.compare("Atom_type:") == 0 ) 		{ ss >> atom_type; }
            if( what.compare("Align:") == 0 )  			{ ss >> mtag_1 >> mtag_2; }

            if( what.compare("Populate:") == 0 )  		{ ss >> population; }

            if( what.compare("Sim_box:") == 0 )             { ss >> sim_box; }

            if( what.compare("Fit") == 0 )  { fit=true; }

            if( what.compare("Impact_vector:") == 0 )  { ss >> ivx.x >> ivx.y >> ivx.z; }

            if( what.compare("Seed:") == 0 )  { ss >> seed; rng.seed(seed); }

            if( what.compare("Lammps_offset:") == 0 )   { ss >> offset; }

            // force-field
            if( what.compare("ff_lj:") == 0 ) { LJ lj; ss >> lj; ff.lj.push_back(lj); }
            if( what.compare("ff_cos2:") == 0 ) { CosSQ cos; ss >> cos; ff.cos.push_back(cos); }

        }
        fs.close();

        return true;
    }

    string toString()
    {
        stringstream ss;

        ss << "Particle_type: " << gen_structure << endl;
        ss << "Output_type: " << out << endl;
        ss << "Number of beads: " << num_of_beads << endl;
        ss << "Scale: " << scale << endl;
        ss << "Offset: " << offset << endl;
        ss << "c: " << c << endl;
        ss << "Number of ligands: " << num_lig << endl;
        ss << "Box: ( " << sim_box << " )" << endl;
        ss << "Position: (" << com_pos.x << ", " << com_pos.y << ", " << com_pos.z << ")" << endl;
        ss << "Patch_1: (" << patch_1.x << "-" << patch_1.vx << ", " << patch_1.y << "-" << patch_1.vy << ", " << patch_1.z << "-" << patch_1.vz << ", " << patch_1.type << ")" << endl;
        ss << "Patch_1: (" << patch_2.x << ", " << patch_2.y << ", " << patch_2.z << ", " << patch_2.type << ")" << endl;
        ss << "Populate: " << population << endl;

        ss << "Force-Field: " << ff << endl;

        return ss.str();
    }

    void clear()
    {
        gen_structure.clear();
        num_of_beads=-1;
        scale=0.0;
        offset=1;
        c=0;
        num_lig=0;
        chain_type=-1;
        mol_tag=-1;
        atom_type=1;

        center=false;
        fit = false;
        mtag_1=-1;
        mtag_2=-1;

        com_pos=Atom(0.0, 0.0, 0.0);

        patches.clear();

        patch_1=Atom(0,0,0,0);
        patch_2=Atom(0,0,0,0);

        ivx=Atom(0.0, 0.0, 0.0);

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

    bool is_mtag_12()
    {
        if( mtag_1 != -1 && mtag_2 != -1 )
            return false;
        return true;
    }

    bool is_mol_tag()
    {
    	return (mol_tag != -1);
    }

    bool isScale()
    {
        return (scale != 0.0);
    }
};


#endif // INPUT_H
