#ifndef DATA_H
#define DATA_H

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
#include <map>

#include "force_field.h"
#include "atom.h"
#include "io_input.h"
#include "io_lammps.h"
#include "io_pdb.h"

#include "virus_pseudot3.h"

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
    vector<Atoms> beads;
    vector<Bonds> bonds;
    vector<Angles> angles;

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
    bool isDefined()
    {
        return !in.infile.empty();
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

            beads.push_back(lammps.beads);
            bonds.push_back(lammps.bonds);
            angles.push_back(lammps.angles);

            lammps.bonds.clear();
            lammps.angles.clear();
            lammps.beads.clear();
        }

        if(in.in.type == IO_Type::pdb)
        {
            pdb.load(in_file);

            beads.push_back(pdb.beads);
            bonds.push_back(Bonds());
            angles.push_back(Angles());

            pdb.beads.clear();
        }

        id_map[in.id] = beads.size()-1;
    }

    ///
    /// Data methods
    ///
    void modify()
    {
        beads[ id_map[in.id] ].scale(in.scale);    // Rescale atom positions
        beads[ id_map[in.id] ].move(in.com_pos);   // Move entire system by vector

        if(in.mol_tag > -1)
            beads[ id_map[in.id] ].set_mol_tag(in.mol_tag); // Change mol_tag of all particles to one set by input

        /*if(in.is_mtag_12())
            align(in.mtag_1, in.mtag_2);*/ // align mol_tag particles in z axis and XY plane

        // if generating into an existing structure that you did not load, give the number of particles as offset
        //offset(all_beads.size());

        if(in.id == 2)
        {
            Virus_pseudoT3 t3;
            t3.set_protomer(beads[0]); // TODO: make it work with input file
            t3.set_capsid(beads[id_map[in.id]]); // TODO: make it work with input file
            t3.calc_protomer_coms();
            t3.calc_pentamers();
            t3.calc_symmetry_axes();

            t3.printXYZ();

            exit(1);
        }

        if( in.center )
            beads[ id_map[in.id] ].center();
    }


    void merge(Atoms& all_beads, Bonds& all_bonds, Angles& all_angles)
    {
        merge(all_beads);
        merge(all_bonds);
        merge(all_angles);
    }

    bool is_overlap(Atoms& other, Force_Field& ff) const
    {
        for(auto& bb : beads)
        {
            if(bb.is_overlap(other, ff))
                return true;
        }
        return false;
    }

    int get_bead_count() const
    {
        int count=0;
        for(auto& bb : beads)
        {
            count += bb.size();
        }
        return count;
    }

    int get_Max_Mol_Tag()
    {
        int max=0;
        for(auto& bb : beads)
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

            lammps.print();
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
        for(auto& item : beads)
        {
            all_beads.insert(all_beads.end(), item.begin(), item.end());
        }
        beads.clear();
    }


    void merge(Bonds& all_bonds)
    {
        for(auto& item : bonds)
        {
            all_bonds.insert(all_bonds.end(), item.begin(), item.end());
        }
        bonds.clear();
    }


    void merge(Angles& all_angles)
    {
        for(auto& item : angles)
        {
            all_angles.insert(all_angles.end(), item.begin(), item.end());
        }
        angles.clear();
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
