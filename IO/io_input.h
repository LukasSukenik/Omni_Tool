#ifndef IO_INPUT_H
#define IO_INPUT_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <cstdlib>
#include <random>
#include <sstream>

#include "sim_box.h"
#include "atom.h"
#include "force_field.h"
#include "rng.h"
#include "types.h"

using namespace std;




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
    unordered_set<string> keys = {"System_type:", "Particle_type:", "System_execute:"};
    unordered_set<string> files = {"Load_file:", "Trajectory_file:", "Trajectory_output_file:", "Histo_2D_dirs_outfile:", "Histo_1D_dirs_outfile:", "Histo_outfile:"};
    Param_Dictionary<string> param = Param_Dictionary<string>(keys + files);

    unordered_set<string> counts1 = {"Number_of_beads:", "Number_of_ligands:", "Num_lipids:", "Number_of_receptors:", "Subdiv_of_beads:", "Subdiv_of_ligands:"};
    unordered_set<string> counts2 = {"Averaged_frame_count:"};
    unordered_set<string> types = {"Mol_tag:", "Chain_type:"};
    unordered_set<string> other = {"ID:", "Seed:"};
    Param_Dictionary<int> p_int = Param_Dictionary<int>(counts1 + counts2 + types + other);

    Param_Dictionary<bool> p_bool = Param_Dictionary<bool>({"Fit:", "Center:", "Only_last_frame:"});
    Param_Dictionary<Tensor_xyz> p_tensor = Param_Dictionary<Tensor_xyz>({"Position_shift:", "Impact_vector:"});
    Param_Dictionary<double> p_float = Param_Dictionary<double>({"Cluster_cutoff:", "Radius:", "Scale:", "b:", "c:", "Cell_size:", "Beads_per_area:", "Ligands_per_area:"});

    Param_Dictionary<vector<int>> p_vec_int = Param_Dictionary<vector<int>>({"Atom_type:", "Atom_mass:", "Histo_settings:", "Histo_spherical_settings:", "Trajectory_settings:"});

    IO out; /// Output type - none, pdb, lammps_full, xyz
    IO in;  /// Input type  - none, pdb, lammps_full
    TMD tmd; // transmembrane domain
    Population population; // Population data
    Simulation_Box sim_box; // simulation box
    vector<LJ> bparam;
    vector<CosSQ> cparam;
    Force_Field ff;

    // Align
    int mtag_1=-1;
    int mtag_2=-1;

    double system_var_a = 0.0;
    double system_var_b = 0.0;
    double bead_size=1.12246204831;

    int offset = 0;

    Atom patch_1 = Atom(0,0,0,0);
    Atom patch_2 = Atom(0,0,0,0);
    vector<Atom> patches;

    // Deprecated - used only within in_input.h
    int op=-1;
    double trim;
    int multiple;

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
        Tensor_xyz v_tensor;
        stringstream ss;

        while( !fs.eof() ) // Lines in input
        {
            key.clear();
            ss.flush();
            ss.clear();
            getline(fs, line);
            ss.str(line);

            ss >> key; // first stuff,

            if(  param.is_key_valid(key) )     { ss >> value;       param[key.substr(0, key.find(':'))] = value; }
            if(  p_bool.is_key_valid(key) )    {                    p_bool[key.substr(0, key.find(':'))] = true; }
            if(  p_int.is_key_valid(key) )     { ss >> value_int;   p_int[key.substr(0, key.find(':'))] = value_int; }
            if(  p_float.is_key_valid(key) )   { ss >> value_float; p_float[key.substr(0, key.find(':'))] = value_float; }
            if(  p_tensor.is_key_valid(key) )  { ss >> v_tensor;    p_tensor[key.substr(0, key.find(':'))] = v_tensor; }
            if(  p_vec_int.is_key_valid(key) ) { p_vec_int[key.substr(0, key.find(':'))] = load_int_array(ss); }

            // Load io
            if( key.compare("Output_type:") == 0 )       { ss >> out; }
            if( key.compare("Input_type:") == 0 )        { ss >> in; }

            // Load particle atom counts
            if( key.compare("Trans_membrane_domain:") == 0 ) { ss >> tmd.size; ss >> tmd.proximal_n; ss >> tmd.distal_n; }

            // Load particle properties: aspect ration, patches
            if( key.compare("Patch:") == 0 )             { Atom a; ss >> a; patches.push_back(a); }

            if( key.compare("Patch_1:") == 0 )           { ss >> patch_1; }
            if( key.compare("Patch_2:") == 0 )           { ss >> patch_2; }

            // Load system
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

            // from gen_membrane - deprecated
            if( key.compare("Operation_type:") == 0 )      { ss >> op; }
            if( key.compare("Trim:") == 0 )                { ss >> trim; }
            if( key.compare("Multiple:") == 0 )            { ss >> multiple; }
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

        for (const auto& [key, value] : param)      { ss << key << ": " << value << '\n'; }
        for (const auto& [key, value] : p_bool)     { ss << key << ": true" << '\n'; }
        for (const auto& [key, value] : p_int)      { ss << key << ": " << value << '\n'; }
        for (const auto& [key, value] : p_float)    { ss << key << ": " << value << '\n'; }
        for (const auto& [key, value] : p_tensor)   { ss << key << ": " << value << '\n'; }
        for (const auto& [key, values] : p_vec_int) { ss << key << ": " << values << '\n'; }

        // Loading a file
        ss << "Input_type: " << in << endl;
        ss << "Output_type: " << out << endl;

        if( param.contains("Particle_type") )
        {
            ss << "Patch_1: (" << patch_1.pos.x << "-" << patch_1.vel.x << ", " << patch_1.pos.y << "-" << patch_1.vel.y << ", " << patch_1.pos.z << "-" << patch_1.vel.z << ", " << patch_1.type << ")" << endl;
            ss << "Patch_2: (" << patch_2.pos.x << ", " << patch_2.pos.y << ", " << patch_2.pos.z << ", " << patch_2.type << ")" << endl;
            ss << "Trans_membrane_domain: (" << tmd.size << ", " << tmd.proximal_n << ", " << tmd.distal_n << ")" << endl;
        }

        if(!sim_box.empty())
            ss << "Box: ( " << sim_box << " )" << endl;

        ss << "Offset: " << offset << endl;

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
        p_bool.clear();
        p_int.clear();
        p_float.clear();
        p_tensor.clear();
        p_vec_int.clear();

        in.clear();
        system_var_a=0.0;
        system_var_b=0.0;

        tmd.clear();

        offset=0;
        mtag_1=-1;
        mtag_2=-1;

        patches.clear();

        patch_1=Atom(0,0,0,0);
        patch_2=Atom(0,0,0,0);

        population.clear();
        bparam.clear();
        cparam.clear();
    }
};


#endif // IO_INPUT_H
