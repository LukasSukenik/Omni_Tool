#ifndef IO_INPUT_H
#define IO_INPUT_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <unordered_map>
#include <variant>

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


enum class IO_Type { none, xyz, pdb, gro, lammps_full };

std::ostream& operator<<(std::ostream& os, const IO_Type out)
{
    switch (out) {
        case IO_Type::gro:
          os << "gro"; return os;
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
        if (str == "gro")
        {
            out.type = IO_Type::gro;
        }
        if (str == "pdb")
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





class TMD
{
public:
    TMD(){}

    int size;
    int proximal_n;
    int distal_n;

    void clear()
    {
        size = 0;
        proximal_n = 0;
        distal_n = 0;
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

    bool empty()
    {
        if(count == 0)
            return true;
        return false;
    }

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


template <typename T> class Param_Dictionary : public unordered_map<string, T>
{
public:
    vector<string> valid_keys;

    Param_Dictionary() {}
};



/**
 * @brief The Input class - Parameters for particle generation
 */
class IO_Input{
public:
    IO_Input() {}

    /**
     * @brief param - List of parameters
     *
     * Load_file: data.start # name of configuration file in format lammps, pdb, or xyz
     * Trajectory_file: traj_1.xtc # name of trajectory file in format xtc, typically file.xtc or traj_[integer].xtc
     * System_type: Flat_Membrane # nane of system: {Flat_Membrane, Lipid_Nanoparticle, Vesicle}
     *
     */
    vector<string> param_valid_key_list = {"Load_file:", "Trajectory_file:", "System_type:"};
    unordered_map<string, string> param;

    vector<string> param_int_valid_key_list = {"ID:"};
    unordered_map<string, int> param_int;

    vector<string> param_float_valid_key_list = {"Cluster_cutoff:"};
    unordered_map<string, double> param_float;

    vector<string> param_vector_int_valid_key_list = {"Atom_type:"};
    unordered_map<string, vector<int>> param_vector_int;

    // IO
    IO out; /// Output type - none, pdb, lammps_full, xyz
    IO in;  /// Input type  - none, pdb, lammps_full

    // particle identifier
    string gen_structure_ID; /// keyword identifying the structure class

    // system identifier
    string system_function;
    double system_var_a = 0.0;
    double system_var_b = 0.0;

    // bead counts
    int num_of_beads=-1; // chain, dodecahedron, ellipsoid, globular_sphere, icosahedron, oblatespheroid, pentamer, slab, sphere, spherepatch, tennisball
    int num_lig=-1;
    double beads_per_area=0.0;
    double ligs_per_area=0.0;
    int subdiv_beads=-1;
    int subdiv_lig=-1;
    TMD tmd; // transmembrane domain

    // from gen_membrane
    int op=-1;
    double radius;
    double trim;
    int num_lipids;
    int multiple;
    int num_rec;

    // atom types
    int chain_type=-1;

    // atom type mass
    vector<int> atom_mass;

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

    // Pore calculation - cell size
    double cell_size=0.0;
    double bead_size=1.12246204831;

    // Population data
    Population population;

    // simulation box
    Simulation_Box sim_box;

    // Force-Field
    vector<LJ> bparam;
    vector<CosSQ> cparam;
    Force_Field ff;

    // Deprecated
    int offset = 0;

    bool is_key_valid(string key)
    {
        for(string valid_key : param_valid_key_list)
        {
            if(key.compare(valid_key) == 0)
            {
                return true;
            }
        }
        return false;
    }

    bool is_key_valid_int(string key)
    {
        for(string valid_key : param_int_valid_key_list)
        {
            if(key.compare(valid_key) == 0)
            {
                return true;
            }
        }
        return false;
    }

    bool is_key_valid_float(string key)
    {
        for(string valid_key : param_float_valid_key_list)
        {
            if(key.compare(valid_key) == 0)
            {
                return true;
            }
        }
        return false;
    }

    bool is_key_valid_vector_int(string key)
    {
        for(string valid_key : param_vector_int_valid_key_list)
        {
            if(key.compare(valid_key) == 0)
            {
                return true;
            }
        }
        return false;
    }


    bool loadInput(string input)
    {
        std::fstream fs( input, std::fstream::in );
        if(!fs.is_open())
        {
            cerr << "Failed to open file: " << input << endl;
            return false;
        }

        string line, key, value;
        int value_int;
        double value_float;
        stringstream ss;
        int len=0;

        while( !fs.eof() ) // Lines in input
        {
            key.clear();
            ss.flush();
            ss.clear();
            getline(fs, line);
            ss.str(line);

            ss >> key; // first stuff,
            is_key_valid(key);

            if(  is_key_valid(key) )            { ss >> value; param[key.substr(0, key.find(':'))] = value; }
            if(  is_key_valid_int(key) )        { ss >> value_int; param_int[key.substr(0, key.find(':'))] = value_int; }
            if(  is_key_valid_float(key) )      { ss >> value_float; param_float[key.substr(0, key.find(':'))] = value_float; }
            if(  is_key_valid_vector_int(key) ) { param_vector_int[key.substr(0, key.find(':'))] = load_atom_type(ss); }

            if( key.compare("Particle_type:") == 0 )     { ss >> gen_structure_ID; }

            // Mixed Params
            if( key.compare("System_execute:") == 0 )    { ss >> system_function >> system_var_a >> system_var_b; }

            // Load io
            if( key.compare("Output_type:") == 0 )       { ss >> out; }
            if( key.compare("Input_type:") == 0 )        { ss >> in; }

            // Load particle atom counts
            if( key.compare("Beads_per_area:") == 0 )   { ss >> beads_per_area; }
            if( key.compare("Ligands_per_area:") == 0 ) { ss >> ligs_per_area; }
            if( key.compare("Number_of_beads:") == 0 )   { ss >> num_of_beads; }
            if( key.compare("Number_of_ligands:") == 0 ) { ss >> num_lig; }
            if( key.compare("Subdiv_of_beads:") == 0 )   { ss >> subdiv_beads; }
            if( key.compare("Subdiv_of_ligands:") == 0 ) { ss >> subdiv_lig; }
            if( key.compare("Trans_membrane_domain:") == 0 ) { ss >> tmd.size; ss >> tmd.proximal_n; ss >> tmd.distal_n; }

            // from gen_membrane
            if( key.compare("Operation_type:") == 0 )      { ss >> op; }
            if( key.compare("Radius:") == 0 )              { ss >> radius; }
            if( key.compare("Trim:") == 0 )                { ss >> trim; }
            if( key.compare("Num_lipids:") == 0 )          { ss >> num_lipids; }
            if( key.compare("Multiple:") == 0 )            { ss >> multiple; }
            if( key.compare("Number_of_receptors:") == 0 ) { ss >> num_rec; }

            // Load atom and molecule types

            if( key.compare("Atom_mass:") == 0 )         { load_atom_mass(ss); }
            if( key.compare("Mol_tag:") == 0 ) 		     { ss >> mol_tag; }
            if( key.compare("Chain_type:") == 0 ) 		 { ss >> chain_type; }

            // Load particle properties: aspect ration, patches
            if( key.compare("b:") == 0 )                 { ss >> b; }
            if( key.compare("c:") == 0 )                 { ss >> c; }
            if( key.compare("Patch:") == 0 )             { Atom a; ss >> a; patches.push_back(a); }
            if( key.compare("Patch_1:") == 0 )           { ss >> patch_1; }
            if( key.compare("Patch_2:") == 0 )           { ss >> patch_2; }

            // Load system
            if( key.compare("Scale:") == 0 )             { ss >> scale; }
            if( key.compare("Center") == 0 )             { center=true; }
            if( key.compare("Position_shift:") == 0 )    { ss >> com_pos.x >> com_pos.y >> com_pos.z; }
            if( key.compare("Seed:") == 0 )              { ss >> seed; rng.seed(seed); }
            if( key.compare("Lammps_offset:") == 0 )     { ss >> offset; }

            // Load calc pore
            if( key.compare("Cell_size:") == 0 )             { ss >> cell_size; }

            // Load the simulation box
            if( key.compare("Sim_box:") == 0 )           { ss >> sim_box; }

            // Load the population of the particle
            if( key.compare("Populate:") == 0 )  		  { ss >> population; }

            // Load the force-field
            if( key.compare("ff_lj:") == 0 )             { LJ lj; ss >> lj; ff.lj[lj.type]=lj; }
            if( key.compare("ff_cos2:") == 0 )           { CosSQ cos; ss >> cos; ff.cos[cos.type]=cos; }

            // Load stuff for nanoparticle orientation and position
            if( key.compare("Align:") == 0 )  			{ ss >> mtag_1 >> mtag_2; }
            if( key.compare("Fit") == 0 )              { fit=true; }
            if( key.compare("Impact_vector:") == 0 )   { ss >> ivx.x >> ivx.y >> ivx.z; }
        }
        fs.close();

        return true;
    }

    vector<int> load_atom_type(stringstream& ss)
    {
        vector<int> atom_type;
        int temp_type;
        while( ss >> temp_type )
        {
            atom_type.push_back(temp_type);
        }
        return atom_type;
    }

    void load_atom_mass(stringstream& ss)
    {
        int temp_mass;
        while( ss >> temp_mass )
        {
            atom_mass.push_back(temp_mass);
        }
    }

    string toString()
    {
        stringstream ss;

        for (const auto& [key, value] : param) {
            ss << key << ": " << value << '\n';
        }

        // Loading a file
        ss << "Input_type:" << in << endl;

        // Generating structure
        if( !gen_structure_ID.empty() )
            ss << "Particle_type: " << gen_structure_ID << endl;

        ss << "Output_type: " << out << endl;

        if( !gen_structure_ID.empty() )
        {
            ss << "Beads_per_area: " << beads_per_area << endl;
            ss << "Ligands_per_area: " << ligs_per_area << endl;
            ss << "Number of beads: " << num_of_beads << endl;
            ss << "Number of ligands: " << num_lig << endl;
            ss << "Subdiv of beads: " << subdiv_beads << endl;
            ss << "Subdiv of ligands: " << subdiv_lig << endl;
            ss << "c: " << c << endl;
            ss << "Patch_1: (" << patch_1.pos.x << "-" << patch_1.vel.x << ", " << patch_1.pos.y << "-" << patch_1.vel.y << ", " << patch_1.pos.z << "-" << patch_1.vel.z << ", " << patch_1.type << ")" << endl;
            ss << "Patch_2: (" << patch_2.pos.x << ", " << patch_2.pos.y << ", " << patch_2.pos.z << ", " << patch_2.type << ")" << endl;
            ss << "Trans_membrane_domain: (" << tmd.size << ", " << tmd.proximal_n << ", " << tmd.distal_n << ")" << endl;
        }

        if(!sim_box.empty())
            ss << "Box: ( " << sim_box << " )" << endl;

        ss << "Scale: " << scale << endl;
        ss << "Offset: " << offset << endl;
        ss << "Position: (" << com_pos.x << ", " << com_pos.y << ", " << com_pos.z << ")" << endl;

        if(!population.empty())
            ss << "Populate: " << population << endl;

        if(!ff.empty())
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

        param.clear();
        param_int.clear();
        param_float.clear();
        param_vector_int.clear();

        gen_structure_ID.clear();
        in.clear();
        system_function.clear();
        system_var_a=0.0;
        system_var_b=0.0;

        beads_per_area=0.0;
        ligs_per_area=0.0;
        num_of_beads=-1;
        num_lig=-1;
        subdiv_beads=-1;
        subdiv_lig=-1;
        tmd.clear();

        scale=1.0;
        offset=0;
        b=0;
        c=0;

        atom_mass.clear();
        chain_type=-1;
        mol_tag=-1;

        center=false;
        fit = false;
        mtag_1=-1;
        mtag_2=-1;

        cell_size = 0.0;

        com_pos=Tensor_xyz(0.0, 0.0, 0.0);

        patches.clear();

        patch_1=Atom(0,0,0,0);
        patch_2=Atom(0,0,0,0);

        ivx=Tensor_xyz(0.0, 0.0, 0.0);

        population.clear();
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
