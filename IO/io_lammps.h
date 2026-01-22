#ifndef IO_LAMMPS_H
#define IO_LAMMPS_H

#include <algorithm>
#include <fstream>
#include <string.h>
#include <array>

#include "types.h"
#include "sim_box.h"




/**
 * @brief The IO_Lammps class
 * Class for IO of files in lammps data formats
 */
class IO_Lammps
{
public:
    IO_Lammps() {}

    bool first_file=true;

    Atoms beads;
    Bonds bonds;
    Angles angles;

    Simulation_Box sim_box;

    vector<My_string > file_head;



    /**
     * Load Lammps full data format file, ignores velocities
     * - in_file - Lammps file name
     */
    void load(string in_file)
    {
        if(in_file.empty()) {
            return;
        }

        cerr << "Loading file " << in_file << endl;
        loadFileHeadAndPart(in_file);
        cerr << "Loading box parameters" << endl;
        loadBox();

        // THROW AWAY VELOCITIES

        cerr << "Parsing bonds parameters" << endl;
        loadBonds(in_file);

        // TEST FOR CONSECCUTIVE INDEXES
        cerr << "Sorting lipids by index" << endl;
        std::sort(beads.begin(), beads.end(), sortN);
        cerr << "Testing index duplicity" << endl;
        for(int i=0; i<beads.size(); i++) {
            if(i%1000 == 0)
            {
                cerr << i << " : " << beads.size() << endl;
            }
            if(i+1 != beads[i].N) {
                cerr << "ERROR, missing index, " << i << " != " << beads[i].N << endl;
            }
        }
        cerr << "Load done" << endl;

        loadAngles(in_file);

    }

    void print()
    {
        cerr << endl;
        cerr << "IO_Lammps::print(), beads:" << beads.size() << ", bonds:" << bonds.size() << ", angles:" << angles.size() << endl;
        if(beads.empty()) {
            cerr << "beads.empty() -> Nothing to generate" << endl;
            return;
        }

        //
        // print force-field
        //
        /*vector<double> dist;
        vector<string> dist_coeff;
        dist = createBondGroups(dist_coeff);
        printForceField(dist, dist_coeff, in.scale);*/

        int bond_types = bonds.calc_Bond_Types();
        int num_a_types = beads.get_Atom_Types().size();

        //
        // Print Head
        //
        cout << "LAMMPS data file via Omni_Tool\n" << endl;
        cout << beads.size() << " atoms\n";
        cout << num_a_types << " atom types\n";

        //
        // Print Bond types
        //
        if(!bonds.empty()) {
            cout << bonds.size() << " bonds\n";
            cout << bond_types << " bond types\n";
        }
        //
        // Print Angle types
        //
        if(!angles.empty()) {
            cout << angles.size() << " angles\n";
            cout << "1 angle types\n";
        }

        //
        // Print Box
        //
        cout << "\n";
        cout << sim_box.xlo << " " << sim_box.xhi << " xlo xhi\n";
        cout << sim_box.ylo << " " << sim_box.yhi << " ylo yhi\n";
        cout << sim_box.zlo << " " << sim_box.zhi << " zlo zhi\n";

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
        for(auto& a : beads) {
            cout << a.N << " " << a.mol_tag << " " << a.type << " " << 0 << " " << a.pos << " 0 0 0" << "\n";
        }

        //
        // Print Bonds
        //
        if(!bonds.empty())
        {
            cout << "\nBonds\n\n";

            for(auto& bond : bonds)
            {
                cout << bond.N << " " << bond.type << " " << bond.at1 << " " << bond.at2 << "\n";
            }
        }

        //
        // Print Angles
        //
        if(!angles.empty()) {
            cout << "\nAngles\n\n";

            for(auto & a : angles) {
                cout << a.N << " " << a.type << " " << a.at1 << " " << a.at2 << " " << a.at3 << "\n";
            }
        }
    }

private:
    void loadFileHeadAndPart(string filename);
    void loadBox();
    void loadBonds(string filename);
    void loadAngles(string filename);

    /*void printForceField(vector<double>& dist, vector<string> &dist_coeff, double scale ) const
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
    }*/

    /*vector<double> createBondGroups(vector<string> &bond_coeff, double precision=1000.0) const
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
    }*/
};




void IO_Lammps::loadAngles(string filename)
{
    char str[256];
    Angle angle;
    std::fstream in;

    in.open(filename, std::fstream::in);
    if(in.good())
    {
        while(in.good())
        {
            in.getline(str, 256);
            if(strstr(str, "Angles") != NULL)
                break;
        }

        while(in.good())
        {
            in >> angle.N >> angle.type >> angle.at1 >> angle.at2 >> angle.at3;
            if( angles.empty() || angle.N != angles.back().N)
            {
                angles.push_back(angle);
            }
        }
    }
    in.close();
}




void IO_Lammps::loadBox()
{
    stringstream str;
    stringstream str2;
    stringstream str3;

    for(unsigned int i=0; i<file_head.size(); i++) {
        //cout << file_head[i].str << endl;
        if(strstr(file_head[i].str, "xlo xhi") != nullptr) {
            str << file_head[i].str;
            str >> sim_box.xlo >> sim_box.xhi;
            continue;
        }

        if(strstr(file_head[i].str, "ylo yhi") != nullptr) {
            str2 << file_head[i].str;
            str2 >> sim_box.ylo >> sim_box.yhi;
            continue;
        }

        if(strstr(file_head[i].str, "zlo zhi") != nullptr) {
            str3 << file_head[i].str;
            str3 >> sim_box.zlo >> sim_box.zhi;
            continue;
        }
    }
}

void IO_Lammps::loadBonds(string filename)
{
    char str[256];
    Bond bond;
    bond.typelock = true;
    std::fstream in;

    in.open(filename, std::fstream::in);
    if(in.good())
    {
        while(!in.eof())
        {
            in.getline(str, 256);
            if(strstr(str, "Bonds") != NULL)
                break;
        }

        cerr << "Found keyword::Bonds" << endl;
        while(!in.eof())
        {
            in.getline(str, 256);
            if (strstr(str, "Angles") != NULL)
                break;

            istringstream iss(str);
            bond = Bond();
            iss >> bond.N >> bond.type >> bond.at1 >> bond.at2;

            if(bond.N > 0)
                bonds.push_back(bond);
        }
    }
    in.close();
    cerr << "File " << filename << " Loaded " << bonds.size() << " bonds" << endl;

    std::sort(bonds.begin(), bonds.end(), Bond::sort_Bond_by_type);

    for(int i=0; i<bonds.size(); i++)
    {
        bonds[i].N = i+1;
        //cerr << temp_bonds[i].N << " " << temp_bonds[i].type << " " << temp_bonds[i].at1 << " " << temp_bonds[i].at2 << endl;
    }
}

void IO_Lammps::loadFileHeadAndPart(string filename)
{
    char str[256];
    Atom part;
    std::fstream in;
    std::stringstream ss;
    int size=beads.size();

    in.open(filename, std::fstream::in);

    if(in.is_open())
    {
        // Load until "Atoms" occurs
        while(in.good() && strstr(str, "Atoms") == NULL)
        {
            in.getline(str, 256);
            if(first_file)
            {
                file_head.push_back(My_string(str));
            }
        }
        in.getline(str, 256); // 1 empty line after atoms

        // Load: N molecule-tag atom-type q x y z nx ny nz
        while(in.good())
        {
            part.N = -1;
            part.nx = 0;
            part.ny = 0;
            part.nz = 0;
            in.getline(str, 256);
            ss.str( str );
            ss >> part.N >> part.mol_tag >> part.type >> part.q >> part.pos.x >> part.pos.y >> part.pos.z >> part.nx >> part.ny >> part.nz;
            ss.flush();
            ss.clear();

            if( part.N == -1 ) {
                if(beads.empty())
                    continue;
                else
                    break;
            }
            beads.push_back(part);
        }
        cerr << "  Added " << beads.size() - size << " beads" << endl;
    } else {
        cerr << "File " << filename << " file not opened" << endl;
    }

    first_file = false;
    in.close();
}

#endif // IO_LAMMPS_H
