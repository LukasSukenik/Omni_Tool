#ifndef IO_INPUT_H
#define IO_INPUT_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <sstream>

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
    unordered_set<string> valid_keys;

    Param_Dictionary(unordered_set<string> valid_key_list) : valid_keys(valid_key_list) {}

    bool is_key_valid(string key) { return valid_keys.contains(key); }
    void validate_keyword(string keyword, string default_value)
    {
        if( !this->contains(keyword) )
        {
            cerr << "Missing keyword; " << keyword << ": " << default_value << endl;
            exit(-1);
        }
    }
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
     * Particle_type: particle identifier::keyword identifying the structure class
     *
     */
    Param_Dictionary<string> param = Param_Dictionary<string>({"Load_file:", "Trajectory_file:", "System_type:", "Particle_type:", "Histo_2D_dirs_outfile:", "Histo_1D_dirs_outfile:"});
    Param_Dictionary<int> p_int = Param_Dictionary<int>({"ID:", "Mol_tag:", "Num_lipids:", "Number_of_receptors:", "Chain_type:", "Seed:", "Averaged_frame_count:"});
    Param_Dictionary<double> p_float = Param_Dictionary<double>({"Cluster_cutoff:", "Radius:", "Scale:", "c:", "Cell_size:"});
    Param_Dictionary<vector<int>> p_vec_int = Param_Dictionary<vector<int>>({"Atom_type:", "Atom_mass:", "Histo_2D_settings:"});

    // IO
    IO out; /// Output type - none, pdb, lammps_full, xyz
    IO in;  /// Input type  - none, pdb, lammps_full

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
    double trim;
    int multiple;

    // molecule types
    int mtag_1=-1;
    int mtag_2=-1;

    // patches
    vector<Atom> patches;
    Atom patch_1 = Atom(0,0,0,0);
    Atom patch_2 = Atom(0,0,0,0);

    // Ellipsoid, oblate spheroid
    myFloat b=0.0;

    // System data
    bool center=false;
    Tensor_xyz com_pos = Tensor_xyz(0.0, 0.0, 0.0);
    Tensor_xyz ivx = Tensor_xyz(0.0, 0.0, 0.0);
    bool fit = false;

    // Pore calculation - cell size
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

        while( !fs.eof() ) // Lines in input
        {
            key.clear();
            ss.flush();
            ss.clear();
            getline(fs, line);
            ss.str(line);

            ss >> key; // first stuff,

            if(  param.is_key_valid(key) )            { ss >> value;       param[key.substr(0, key.find(':'))] = value; }
            if(  p_int.is_key_valid(key) )        { ss >> value_int;   p_int[key.substr(0, key.find(':'))] = value_int; }
            if(  p_float.is_key_valid(key) )      { ss >> value_float; p_float[key.substr(0, key.find(':'))] = value_float; }
            if(  p_vec_int.is_key_valid(key) ) { p_vec_int[key.substr(0, key.find(':'))] = load_int_array(ss); }

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
            if( key.compare("Trim:") == 0 )                { ss >> trim; }
            if( key.compare("Multiple:") == 0 )            { ss >> multiple; }

            // Load particle properties: aspect ration, patches
            if( key.compare("b:") == 0 )                 { ss >> b; }
            if( key.compare("Patch:") == 0 )             { Atom a; ss >> a; patches.push_back(a); }
            if( key.compare("Patch_1:") == 0 )           { ss >> patch_1; }
            if( key.compare("Patch_2:") == 0 )           { ss >> patch_2; }

            // Load system
            if( key.compare("Center") == 0 )             { center=true; }
            if( key.compare("Position_shift:") == 0 )    { ss >> com_pos.x >> com_pos.y >> com_pos.z; }
            if( key.compare("Lammps_offset:") == 0 )     { ss >> offset; }



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

        if(p_int.contains("Seed"))
        {
            rng.seed(p_int["Seed"]);
        }

        return true;
    }

    vector<int> load_int_array(stringstream& ss)
    {
        vector<int> int_array;
        int temp;
        while( ss >> temp )
        {
            int_array.push_back(temp);
        }
        return int_array;
    }

    string toString()
    {
        stringstream ss;

        for (const auto& [key, value] : param) {
            ss << key << ": " << value << '\n';
        }
        for (const auto& [key, value] : p_int) {
            ss << key << ": " << value << '\n';
        }
        for (const auto& [key, value] : p_float) {
            ss << key << ": " << value << '\n';
        }
        for (const auto& [key, values] : p_vec_int) {
            ss << key << ": ";
            for(const auto& val : values)
            {
                ss << val << ' ';
            }
            ss << '\n';
        }


        // Loading a file
        ss << "Input_type:" << in << endl;

        ss << "Output_type: " << out << endl;

        if( param.contains("Particle_type") )
        {
            ss << "Beads_per_area: " << beads_per_area << endl;
            ss << "Ligands_per_area: " << ligs_per_area << endl;
            ss << "Number of beads: " << num_of_beads << endl;
            ss << "Number of ligands: " << num_lig << endl;
            ss << "Subdiv of beads: " << subdiv_beads << endl;
            ss << "Subdiv of ligands: " << subdiv_lig << endl;
            ss << "Patch_1: (" << patch_1.pos.x << "-" << patch_1.vel.x << ", " << patch_1.pos.y << "-" << patch_1.vel.y << ", " << patch_1.pos.z << "-" << patch_1.vel.z << ", " << patch_1.type << ")" << endl;
            ss << "Patch_2: (" << patch_2.pos.x << ", " << patch_2.pos.y << ", " << patch_2.pos.z << ", " << patch_2.type << ")" << endl;
            ss << "Trans_membrane_domain: (" << tmd.size << ", " << tmd.proximal_n << ", " << tmd.distal_n << ")" << endl;
        }

        if(!sim_box.empty())
            ss << "Box: ( " << sim_box << " )" << endl;

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
        p_int.clear();
        p_float.clear();
        p_vec_int.clear();

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

        offset=0;
        b=0;

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
        bparam.clear();
        cparam.clear();
    }
};


#endif // IO_INPUT_H
