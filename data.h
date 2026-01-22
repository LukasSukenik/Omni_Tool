#ifndef DATA_H
#define DATA_H

#include <iostream>
#include <fstream>
#include <sstream>

#include <iomanip>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <random>
#include <string.h>

#include <array>
#include <map>

#include "force_field.h"
#include "atom.h"
#include "io_input.h"
#include "io_lammps.h"
#include "io_pdb.h"

//#include "virus_pseudot3.h"

using namespace std;




/**
 * @brief The Data class - Generated data from class particle
 */
class Data
{
public:
	//
    /// Data
	//
    map<int, int> id_map;
    vector<Atoms> coll_beads;
    vector<Bonds> coll_bonds;
    vector<Angles> coll_angles;

    //
    /// force field stuff - deprecated
    //
    array<array<double, 100>, 100> all_sigma;
    int all_sigma_size = 0;
    array< array<bool, 100>, 100> all_sigma_cosatt;
    vector<LJ> all_bparam;
    vector<CosSQ> all_cparam;

    //
    // Class for loading input file
    //
    IO_Input in;
    IO_Lammps lammps;
    IO_PDB pdb;
    IO_XYZ xyz;



    Data()
    {
        for(int i=0; i<99; ++i)
            for(int j=0; j<99; ++j)
                all_sigma_cosatt[i][j] = false;
        for(int i=0; i<99; ++i)
            for(int j=0; j<99; ++j)
                all_sigma[i][j] = 0.0;
    }

    ///
    /// Input methods
    ///
    bool is_load_file()
    {
        return !in.file_structure.empty();
    }

    bool is_particle_gen()
    {
        return !in.gen_structure_ID.empty();
    }

    bool is_system()
    {
        return !in.system_type.empty();
    }

    bool load_input(string input)
    {
        in.clear();
        in.loadInput(input);
        return true;
    }

    void load_data(string in_file)
    {
        if(in.in.type == IO_Type::lammps_full)
        {
            lammps.load(in_file);
            in.sim_box = lammps.sim_box;

            coll_beads.push_back(lammps.beads);
            coll_bonds.push_back(lammps.bonds);
            coll_angles.push_back(lammps.angles);

            lammps.bonds.clear();
            lammps.angles.clear();
            lammps.beads.clear();
        }

        if(in.in.type == IO_Type::pdb)
        {
            pdb.load(in_file);

            coll_beads.push_back(pdb.beads);
            coll_bonds.push_back(Bonds());
            coll_angles.push_back(Angles());

            pdb.beads.clear();
        }

        id_map[in.id] = coll_beads.size()-1;

        if(coll_beads.empty())
        {
            cerr << "Data::load_data coll_beads is empty, nothing was loded" << endl;
            exit(1);
        }
    }

    ///
    /// Data methods
    ///
    void modify()
    {
        coll_beads[ id_map[in.id] ].scale(in.scale);    // Rescale atom positions
        coll_beads[ id_map[in.id] ].move(in.com_pos);   // Move entire system by vector

        if(in.mol_tag > -1)
            coll_beads[ id_map[in.id] ].set_mol_tag(in.mol_tag); // Change mol_tag of all particles to one set by input

        /*if(in.is_mtag_12())
            align(in.mtag_1, in.mtag_2); // align mol_tag particles in z axis and XY plane
        */

        // if generating into an existing structure that you did not load, give the number of particles as offset
        //offset(all_beads.size());

        if( in.center )
            coll_beads[ id_map[in.id] ].center();
    }


    void merge(Atoms& all_beads, Bonds& all_bonds, Angles& all_angles)
    {
        merge(all_beads);
        merge(all_bonds);
        merge(all_angles);
    }

    bool is_overlap(Atoms& other, Force_Field& ff) const
    {
        for(auto& bb : coll_beads)
        {
            if(bb.is_overlap(other, ff))
                return true;
        }
        return false;
    }

    int get_bead_count() const
    {
        int count=0;
        for(auto& bb : coll_beads)
        {
            count += bb.size();
        }
        return count;
    }

    int get_bond_count() const
    {
        int count=0;
        for(auto& bb : coll_bonds)
        {
            count += bb.size();
        }
        return count;
    }

    int get_Max_Mol_Tag()
    {
        int max=0;
        for(auto& bb : coll_beads)
        {
            if(max < bb.get_Max_Mol_Tag())
            {
                max = bb.get_Max_Mol_Tag();
            }
        }
        return max;
    }

    ///
    /// Output methods
    ///
    void print()
    {
        Atoms all_beads;
        Bonds all_bonds;
        Angles all_angles;

        merge(all_beads, all_bonds, all_angles);

        if( in.out.type == IO_Type::lammps_full)
        {
            lammps.sim_box = in.sim_box;

            lammps.beads = all_beads;
            lammps.bonds = all_bonds;
            lammps.angles = all_angles;

            lammps.print(in);
        }
        if( in.out.type == IO_Type::xyz)
        {
            xyz.beads = all_beads;
            xyz.print();
        }
        if( in.out.type == IO_Type::pdb)
        {
            pdb.beads = all_beads;
            pdb.print();
        }

        report(all_beads);
    }

    void report(Atoms all_beads)
    {
        if(all_beads.empty())
        {
            cerr << "!!! Nothing generated !!!" << endl;
        }
        else
        {
            cerr << all_beads.toString() << endl;
        }
    }

private:
    void merge(Atoms& all_beads)
    {
        for(auto& item : coll_beads)
        {
            all_beads.insert(all_beads.end(), item.begin(), item.end());
        }
        coll_beads.clear();
    }


    void merge(Bonds& all_bonds)
    {
        for(auto& item : coll_bonds)
        {
            all_bonds.insert(all_bonds.end(), item.begin(), item.end());
        }
        coll_bonds.clear();
    }


    void merge(Angles& all_angles)
    {
        for(auto& item : coll_angles)
        {
            all_angles.insert(all_angles.end(), item.begin(), item.end());
        }
        coll_angles.clear();
    }

    LJ getBeadParam(int type)
    {
        for(auto item : all_bparam)
        {
            if( type == item.type)
            {
                return item;
            }
        }
        return LJ();
    }
};




#endif // DATA_H
