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

#include "force_field.h"
#include "atom.h"
#include "io_input.h"
#include "io_lammps.h"
#include "io_pdb.h"

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

    bool load_input(string input, int i)
    {
        in.clear();
        in.loadInput(input, i);
        return true;
    }

    void load_data(string in_file)
    {
        if(in.in.type == IO_Type::lammps_full)
        {
            lammps.load(in_file);
            in.sim_box = lammps.sim_box;
            bonds.push_back(lammps.bonds);
            angles.push_back(lammps.angles);
            beads.push_back(lammps.beads);

            lammps.bonds.clear();
            lammps.angles.clear();
            lammps.beads.clear();
        }

        if(in.in.type == IO_Type::pdb)
        {
            pdb.load(in_file);
            beads.push_back(pdb.beads);
            pdb.beads.clear();
        }
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

    ///
    /// Data methods
    ///
    void modify()
    {
        cout << in.i << endl;
        exit(-1);

        beads[in.i].scale(in.scale);    // Rescale atom positions
        beads[in.i].move(in.com_pos);   // Move entire system by vector

        if(in.mol_tag > -1)
            beads[in.i].set_mol_tag(in.mol_tag); // Change mol_tag of all particles to one set by input
        /*if(in.is_mtag_12())
            align(in.mtag_1, in.mtag_2);*/ // align mol_tag particles in z axis and XY plane

        // if generating into an existing structure that you did not load, give the number of particles as offset
        //offset(all_beads.size());


        if(beads.size() == 2)
        {
            Atoms sym_axis_2;
            Atoms sym_axis_3;
            Atoms sym_axis_5;

            identify_symmetry(sym_axis_2, sym_axis_3, sym_axis_5);

            align_2fold_on_x();
        }

        if( in.center )
            beads[in.i].center();
    }

    void merge(Atoms& all_beads, Bonds& all_bonds, Angles& all_angles)
    {
        merge(all_beads);
        merge(all_bonds);
        merge(all_angles);
    }

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

private:
    void align_2fold_on_x()
    {
        if(beads.size() == 2)
        {
            ;
            //cout << "exit in Data::align_2fold_on_x" << endl;
            //exit(1);
        }
    }


    void identify_symmetry(Atoms& sym_axis_2, Atoms& sym_axis_3, Atoms& sym_axis_5)
    {
        Atoms com = get_protomer_coms();
        sym_axis_5 = get_5_axis(com);
        sym_axis_2 = get_2_axis(sym_axis_5);
        sym_axis_3 = get_3_axis(sym_axis_5);

        beads.push_back(sym_axis_2);
        beads.push_back(sym_axis_3);
        beads.push_back(sym_axis_5);

        beads[0].clear();
        //beads[1].clear();
    }


    Atoms get_protomer_coms()
    {
        int proto_size = beads[0].size();
        Atoms com;
        Atom cm;

        // Generate protomers COMs
        for(int i=0; i<60; ++i)
        {
            cm = beads[1].center_of_mass( i*proto_size, (i+1)*proto_size );
            com.push_back(cm);
        }

        return com;
    }

    Atoms get_2_axis(Atoms& sym_axis_5)
    {
        bool exists = false;
        double min_d = sym_axis_5.min_dist();
        Atoms sym_axis_2;

        for(auto& i : sym_axis_5)
        {
            for(auto& j : sym_axis_5)
            {
                if(i != j && i.dist(j) < min_d*1.05)
                {
                    Atom add = i+j;
                    add = add * (1.0/add.size()) * i.size() ;
                    add.atom_name = " N  ";
                    add.element = " N";

                    exists = false;
                    for(auto& a : sym_axis_2)
                    {
                        if(a.isAproxSame(add))
                            exists = true;
                    }

                    if(!exists)
                        sym_axis_2.push_back( add );
                }
            }
        }

        return sym_axis_2;
    }

    Atoms get_3_axis(Atoms& sym_axis_5)
    {
        bool exists = false;
        double min_d = sym_axis_5.min_dist();
        Atoms sym_axis_3;

        for(auto& i : sym_axis_5)
        {
            for(auto& j : sym_axis_5)
            {
                for(auto& k : sym_axis_5)
                {
                    if( i != j && i.dist(j) < min_d*1.05 &&
                        i != k && i.dist(k) < min_d*1.05 &&
                        j != k && j.dist(k) < min_d*1.05 )
                    {
                        Atom add = i+j+k;
                        add = add * (1.0/add.size()) * i.size() ;
                        add.atom_name = " O  ";
                        add.element = " O";

                        exists = false;
                        for(auto& a : sym_axis_3)
                        {
                            if(a.isAproxSame(add))
                                exists = true;
                        }

                        if(!exists)
                            sym_axis_3.push_back( add );
                    }
                }
            }
        }

        return sym_axis_3;
    }

    Atoms get_5_axis(Atoms& com)
    {
        bool exists = false;
        Atoms sym_axis_5;
        Atoms sel;
        vector<Atoms> penta;

        for(auto& item : com)
        {
            sel = get_pentamer(com, item);
            if(penta.empty())
            {
                penta.push_back(sel);
                sym_axis_5.push_back( sel.center_of_mass() );
            }

            exists = false;
            for(auto& p : penta)
            {
                if( sel.similar(p) )
                {
                    exists = true;
                }
            }
            if(!exists)
            {
                penta.push_back(sel);
                sym_axis_5.push_back( sel.center_of_mass() );
            }
        }

        sym_axis_5.scale( com[0].size() / sym_axis_5[0].size() );

        return sym_axis_5;
    }


    Atoms get_pentamer(Atoms& container, Atom origin)
    {
        Atoms sel;
        Atoms penta;
        double min_d = container.min_dist();

        for(double i : {1.5, 1.6, 1.7, 1.8, 1.9})
        {
            sel = within(origin, container, i*min_d);
            if(sel.size() > 5)
                break;
        }

        // Remove extraneous protomer COMS
        vector<int> nei;

        for(auto& item : sel)
        {
            nei.push_back( within(item, sel, min_d*1.01).size() );
        }

        for(int i=0; i<sel.size(); ++i)
        {
            if(nei[i] >= 3)
            {
                penta.push_back(sel[i]);
            }
        }

        return penta;
    }



    Atoms within(Atom& start, Atoms& container, double radius)
    {
        Atoms ret;

        for(auto& a : container)
        {
            if( start.dist(a) < radius )
            {
                ret.push_back(a);
            }
        }

        return ret;
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
