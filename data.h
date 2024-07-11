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
    /// lammps data stuff
	//
	Atoms all_beads;
    Bonds all_bonds;
    Angles all_angles;

    vector<LJ> all_bparam;
    vector<CosSQ> all_cparam;

    Atoms temp_beads;
    Bonds temp_bonds;
    Angles temp_angles;

    //
    /// force field stuff
    //
    array<array<double, 100>, 100> all_sigma;
    int all_sigma_size = 0;
    array< array<bool, 100>, 100> all_sigma_cosatt;





    //
    // Class for loading input file
    //
    Input in;
    IO_Lammps lammps;
    IO_PDB pdb;



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
        lammps.load(in_file);
        in.sim_box = lammps.sim_box;
        temp_bonds = lammps.bonds;
        temp_angles = lammps.angles;
        temp_beads = lammps.beads;
    }

    ///
    /// Output methods
    ///

    void printForceField(vector<double>& dist, vector<string> &dist_coeff, double scale ) const
    {
        fstream force_field("force_field", fstream::out);

        force_field << "variable coeff_1 string 1.0" << endl;
        force_field << "variable coeff_2 string 1.0" << endl;
        force_field << "variable coeff_3 string 1.0" << endl;
        force_field << "variable coeff_bond string 50" << endl;
        force_field << endl;

        for(int i=0; i<all_sigma_size; ++i) {
            for(int j=i; j<all_sigma_size; ++j) {
                force_field << "pair_coeff " << i+1 << " " << j+1 << " lj/cut 1.0 " << 0.5 * (all_sigma[i][j]+all_sigma[i][j]) / 1.122462048 << " " << 0.5 * (all_sigma[i][j]+all_sigma[i][j]) << endl;
            }
            force_field << endl;
        }

        force_field << endl;
        int count = 1;
        for(int i=0; i<all_sigma_size; ++i) {
            for(int j=i; j<all_sigma_size; ++j) {
                if(all_sigma_cosatt[i][j]) {
                    force_field << "pair_coeff " << i+1 << " " << j+1 << " cosatt ${coeff_" << count << "} " << 0.5 * (all_sigma[i][j]+all_sigma[i][j]) << " 1.0"<< endl; // 2 3
                    ++count;
                }
            }
        }

        force_field << "\n" << endl;

        for(int j=0; j<dist.size(); ++j) {
            force_field << "bond_coeff " << j+1 << " harmonic ${" << dist_coeff[j] << "} " << dist[j] << endl;
        }

        force_field.close();
    }

    void print() const
    {
        if( in.out.type == Output_Type::lammps_full)
        {
            printLammps();
        }
        if( in.out.type == Output_Type::xyz)
        {
            printXYZ();
        }
        if( in.out.type == Output_Type::pdb)
        {
            printPDB();
        }
    }

    void printLammps() const
    {
        if(all_beads.empty()) {
            cerr << "No beads generated" << endl;
            return;
        }

        vector<double> dist;
        vector<string> dist_coeff;

        dist = createBondGroups(dist_coeff);
        printForceField(dist, dist_coeff, in.scale);

        int bond_types = all_bonds.calc_Bond_Types();
        int num_a_types = all_beads.get_Atom_Types().size();

        //
        // Print Head
        //
        cout << "LAMMPS data file via Omni_Tool\n" << endl;
        cout << all_beads.size() << " atoms\n";
        cout << num_a_types << " atom types\n";

        //
        // Print Bond types
        //
        if(!all_bonds.empty()) {
            cout << all_bonds.size() << " bonds\n";
            cout << bond_types << " bond types\n";
        }
        //
        // Print Angle types
        //
        if(!all_angles.empty()) {
            cout << all_angles.size() << " angles\n";
            cout << "1 angle types\n";
        }

        //
        // Print Box
        //
        cout << "\n";
        cout << in.sim_box.xlo << " " << in.sim_box.xhi << " xlo xhi\n";
        cout << in.sim_box.ylo << " " << in.sim_box.yhi << " ylo yhi\n";
        cout << in.sim_box.zlo << " " << in.sim_box.zhi << " zlo zhi\n";

        //
        // Print Masses
        //
        cout << "\nMasses\n\n";
        for(int i=0; i<num_a_types; ++i )
            cout << i+1 << " mass_" << i+1 << "\n";

        //
        // Print Atoms
        //
        cout <<"\nAtoms # full\n" << endl;
        for(auto& a : all_beads) {
            cout << a.N << " " << a.mol_tag << " " << a.type << " " << 0 << " " << a.pos << " 0 0 0" << "\n";
        }

        //
        // Print Bonds
        //
        if(!all_bonds.empty())
        {
            cout << "\nBonds\n\n";

            for(auto& bond : all_bonds)
            {
                cout << bond.N << " " << bond.type << " " << bond.at1 << " " << bond.at2 << "\n";
            }
        }

        //
        // Print Angles
        //
        if(!all_angles.empty()) {
            cout << "\nAngles\n\n";

            for(auto & a : all_angles) {
                cout << a.N << " " << a.type << " " << a.at1 << " " << a.at2 << " " << a.at3 << "\n";
            }
        }
    }

    void printXYZ() const
    {
        cout << all_beads.size() << "\nparticle\n";
        for (const Atom& atom : all_beads)
        {
            cout << "C" << atom.type <<  " " << atom.pos << "\n";
        }
    }

    void printPDB() const
    {
        vector<string> type_to_atom_name;
        type_to_atom_name.push_back("H");
        type_to_atom_name.push_back("B");
        type_to_atom_name.push_back("C");
        type_to_atom_name.push_back("N");
        type_to_atom_name.push_back("O");
        type_to_atom_name.push_back("F");
        type_to_atom_name.push_back("P");
        type_to_atom_name.push_back("K");

        for (const Atom& atom : all_beads) {
            cout << "ATOM" << "  " // 1-4, 5-6:empty
                 << std::setw(5) << atom.N << " " // 7-11: atom serial number, 12:empty
                 << std::setw(4) << std::left << type_to_atom_name[atom.type] << " " // 13-16:atom name, 17:empty
                 << std::setw(3) << std::right << atom.type << " " // 18-20: Residue name, 21:empty
                 << std::setw(1) << atom.mol_tag // 22: chain identifier
                 << std::setw(4) << atom.type // 23-26: Residue sequence number
                 << " " << "   " // 27:code for insertion of residues, 28-30:empty
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.pos.x // 31-38
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.pos.y // 39-46
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.pos.z // 47-54
                 << "   1.0" // 55-60: Occupancy
                 << "   1.0" // 61-66: Temperature factor
                 << "      " // 67-72: Empty
                 << "    " // 73-76: Segment identifier (optional)
                 << std::setw(2) << type_to_atom_name[atom.type] // 77-78 Element symbol
                 << "  " << "\n"; // 79-80 Charge (optional)
          }

    }

    ///
    /// Data methods
    ///
    void modify()
    {
        scale(in.scale);    // Rescale atom positions
        move(in.com_pos);   // Move entire system by vector

        if(in.is_mol_tag())
            set_mol_tag(in.mol_tag); // Change mol_tag of all particles to one set by input
        if(in.is_mtag_12())
            align(in.mtag_1, in.mtag_2); // align mol_tag particles in z axis and XY plane

        // if generating into an existing structure that you did not load, give the number of particles as offset
        offset(all_beads.size());

        if( in.fit )
            fit();
        if( in.center )
            center();
    }

    void add()
    {
        all_beads.insert(all_beads.end(), temp_beads.begin(), temp_beads.end());
        all_bonds.insert(all_bonds.end(), temp_bonds.begin(), temp_bonds.end());
        all_angles.insert(all_angles.end(), temp_angles.begin(), temp_angles.end());

        temp_beads.clear();
        temp_bonds.clear();
        temp_angles.clear();
    }




//private:

    void offset(int offs)
    {
    	temp_beads.offset(offs);
    	temp_bonds.offset(offs);
    	temp_angles.offset(offs);
    }

    void move(Atom move)
    {
    	temp_beads.move(move);
    	cerr << "move " << move.pos.x << " " << move.pos.y << " " << move.pos.z << " done" << endl;
    }

    /**
     * @brief center - moves structure, so that center_of_mass is (0,0,0)
     * @param mtag
     */
    void center(int mtag=-1)
    {
        if(! temp_beads.empty())
        {
            Atom cm = temp_beads.center_of_mass(mtag);
            cm*=-1.0;
            move( cm );
            cerr << "center of " << mtag << " done" << endl;
        }
        else
        {
            cerr << "Error: No atoms loaded" << endl;
        }
    }

    void scale(double scale)
    {
    	if(in.scale != 0.0)
    	{
    		temp_beads.scale(scale);
    	}
    	cerr << "scale " << scale << " done" << endl;
    }

    void set_mol_tag(int mtag)
    {
        temp_beads.set_mol_tag(mtag);
        cerr << "All beads changed to mol_tag = " << mtag << endl;
    }

    /**
     * @brief impact - Deprecated, don't use
     * @param ivx
     */
    void impact(Tensor_xyz ivx)
    {
        if(ivx.size() > 0.1)
        {


            // move the liposome to COM
            Atom com = temp_beads.center_of_mass();
            com *= -1;
            move(com);

            // Define nanoparticle and liposomes dimensions
            double z_min_lipo3 = 9999; // liposome we are adding
            double y_min_lipo3 = 9999; // liposome we are adding
            double y_max_lipo3 = -9999; // liposome we are adding
            double z_max_nano = -9999;  // nanoparticle
            double z_max_lipo1 = -9999; // bound liposome
            double y_max_lipo1 = -9999; // bound liposome
            double y_min_lipo1 = 9999; // bound liposome


            for(Atom& item : temp_beads)
            {
                if(z_min_lipo3 > item.pos.z)
                    z_min_lipo3 = item.pos.z;
                if(y_min_lipo3 > item.pos.y)
                    y_min_lipo3 = item.pos.y;
                if(y_max_lipo3 < item.pos.y)
                    y_max_lipo3 = item.pos.y;
            }

            for(Atom& item : all_beads)
            {
                if(y_min_lipo1 > item.pos.y && item.mol_tag == 1)
                    y_min_lipo1 = item.pos.y;
                if(y_max_lipo1 < item.pos.y && item.mol_tag == 1)
                    y_max_lipo1 = item.pos.y;
                if(z_max_lipo1 < item.pos.z && item.mol_tag == 1)
                    z_max_lipo1 = item.pos.z;
                if(z_max_nano < item.pos.z && item.mol_tag == 2)
                    z_max_nano = item.pos.z;
            }

            // z impact
            //z_impact.z = z_max_nano - z_min_lipo3 + 3;

            // y impact
            //y_impact.y = - y_max_lipo3 - z_max_nano -1;
            //y_impact.z = z_max_lipo1 - z_min_lipo3 + 0.5;

            ivx.normalise();

            Atom impact = Atom(0.0, 0.0, 0.0);
            impact.pos.y = (- y_max_lipo3 - z_max_nano -1) * ivx.y;
            impact.pos.z = (z_max_lipo1 - z_min_lipo3 + 0.5)*(1-ivx.z) + ivx.z*(z_max_nano - z_min_lipo3 + 3);
            move(impact);

            cerr << "Liposome moved by " << impact.pos.x << " " << impact.pos.y << " " << impact.pos.z << endl;
        }
    }

    /**
     * @brief fit positions the loaded/generated structure next to previosly generated/loadedd structure
     * - used for ideal collision position of two liposomes and a nanoparticle
     * see Fit_function.blend, need blender 2.9
     */
    void fit()
    {
        //
        // TODO: Construct an algorithm for positioning the second liposome in an ideal collision position.
        // - I think something as shown in Fit_function.blend will work fine, but feel free to innovate
        //
        // Moving the structure is already implemented in move(Atom vec) function, example below
        // Calculating overlap is simple as well, example below
        //

        Atom displace = Atom(0, 0, 23); // class Atom works as a vector as well.
        move(displace); // displace liposome2 by vector displace
        return;

        //for(Atom& item : temp_beads) // loop over liposome2
        //for(Atom& item : all_beads) // loop over liposome+nanoparticle structure
        double distance_squared; // we are using distance squared because it takes less resources -> we are not calculating the square root
        double too_small = 1.4; // maybe bigger, smaller? Best to eyeball it for vmd once you make it semi-functional
        for(Atom& lip2 : temp_beads)
        {
            for(Atom& nano_lip : all_beads)
            {
                if(nano_lip.mol_tag == 2) // 2 is nanoparticle in our use case, 1 is the liposome
                {
                    distance_squared = lip2.distSQ(nano_lip); // assume overlap if distance squared too small
                    if(distance_squared < too_small)
                    {
                        cerr << "Overlap!" << endl;
                    }
                }
            }
        }

        //
        // Rotating a structure is shown in align function
        //

        //
        // Finally Calculate impact vector and print it out into a file, then load in in prep.sh script
        // tutorial for input/output in c++ https://www.cplusplus.com/doc/tutorial/files/
        //
    }

    /**
     * @brief align - align a liposome and nanoparticle
     * @param mtag - nanoparticle mol_tag
     * @param mtag2 - liposome mol_tag
     */
    void align(int mtag, int mtag2)
    {
        // Test for empty mol_tags
        Atom x_axis = Atom(1,0,0);
        Atom x_axis_negative = Atom(-1,0,0);
        Atom z_axis = Atom(0,0,-1);

        //
        // Move nanoparticle to center (0,0,0)
        //
        center(mtag);

        //
        // Rotate structure so that mtag beads (nanoparticle) align with x_axis
        // - nanoparticle generated from poles = tips in prolate form, same as in oblate form
        // -- 1/4 beads from each end identify the poles (tips)
        //
        int count = temp_beads.count_atoms_of_Mol_tag(mtag);             // number of mtag beads(nanoparticle beads)
        Atom nano1 = temp_beads.center_of_mass(mtag, 0, count/4);         // first 1/4 COM of mtag beads
        Atom nano2 = temp_beads.center_of_mass(mtag, 1+3*count/4, count); // last 1/4 COM of mtag beads
        Atom nano_axis = nano1-nano2;                          // Axis of mtag beads
        nano_axis.normalise();                                 // normalise axis vector for correct rotation
        //
        // Look at vector_magic.blend for visual example, need blender 2.9
        // - rot_axis defined plane is between the nano_axis and x_axis vector
        // - by rotating the nano_axis vector 180deg in this plane we align in at x_axis_negative exactly
        //
        Atom rot_axis = nano_axis-x_axis;
        rot_axis.normalise();                                  // normalise for rotation algo
        for(Atom& item : temp_beads)
        {
            item.rotate(rot_axis, 3.14159265359);       // rotate 180deg = 3.1415 radians, rad to deg = 57.2958
        }

        //
        // Rotate mtag so that COM of mtag2 is (*,*,0) = centered around z axis
        // - to keep it aligned with x axis, we rotate only around x axis
        //
        Atom com_mtag2 = temp_beads.center_of_mass(mtag2);
        com_mtag2.pos.x = 0.0;
        com_mtag2.normalise();
        double angle = acos( com_mtag2.dot(z_axis) );
        double clockwise = (com_mtag2.cross(z_axis)).pos.x ;

        if(clockwise > 0.0)
        {
            for(Atom& item : temp_beads)
            {
                item.rotate(x_axis, angle);       //
            }
        } else
        {
            for(Atom& item : temp_beads)
            {
                item.rotate(x_axis_negative, angle);       //
            }
        }

        //
        // TODO: Construct nanoparticle patch, rotate structure so patch points to in +y axis
        //
        /*Atom patch_vec;
            Atom rotate_axis;
            double rotate_angle;*/

        //
        // Class Atom has variable x,y,z. You can access them via . (dot)
        // example: patch_vec.y
        //
        // there are a number of function within class Atom that you can use. -, +, *, /, dot, cross.
        // If you are unsure what they do look at how they are programmed in class Atom
        //

        //
        // function center_of_mass(mol_tag) returns position of center of mass. This is stored in class Atom
        //
        //Atom mtag2_COM = center_of_mass(mtag2);

        //
        // Rotates structure around rotate_axis by angle rotate_angle
        //
        /*rotate_axis.normalise();
            for(Atom& item : temp_beads)
            {
                item.rotate(rotate_axis, rotate_angle);       //
            }*/

        cerr << "Aligned to x axis and z axis" << endl;
    }



    //
    // all_beads functions
    //

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

    string toString()
    {
        stringstream ss;
        vector<int> moltags = getMolTypes();
        vector<int> types = all_beads.get_Atom_Types();

        ss << "Beads: " << all_beads.size() << endl;
        ss << types.size() << " atom types:" << endl;
        for(int atom_type : types)
        {
            ss << "Atom type " << atom_type << " of " << all_beads.count_Atoms_of_Type(atom_type) << endl;
        }

        ss << moltags.size() << " molTypes:" << endl;
        for(int mol_tag : moltags)
        {
            ss << "Molecule " << mol_tag << " with " << all_beads.count_atoms_of_Mol_tag(mol_tag) << " atoms" << endl;
        }

        return ss.str();
    }

    vector<int> getMolTypes()
    {
    	return all_beads.get_Mol_Types();
    }



    void roundBondDist(vector<double>& dist, double precision=1000.0)
    {
        for(int j=0; j<dist.size(); ++j) // each j+1 is bond type
        {
            for(auto& bond : all_bonds)
            {
                if( !bond.typelock && isAproxSame( bond.r0, dist[j], 1.0/precision) )
                {
                    bond.type = j+1;
                }
            }
        }
    }



    vector<double> createBondGroups(vector<string> &bond_coeff, double precision=1000.0) const
    {
        vector<double> dist;
        if(!all_bonds.empty())
        {
            bool exist;
            for(auto& bond : all_bonds)
            {
                if( !bond.typelock )
                {
                    exist = false;
                    for(int j=0; j<dist.size(); ++j)
                    {
                        if( isAproxSame( bond.r0, dist[j], 1.0/precision) )
                        {
                            exist=true;
                        }
                    }
                    if(!exist)
                    {
                        dist.push_back( round(precision * bond.r0) / precision );
                        bond_coeff.push_back( bond.coeff_name );
                    }
                }
            }
        }
        return dist;
    }


};




#endif // DATA_H
