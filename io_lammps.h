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

private:
    void loadFileHeadAndPart(string filename);
    void loadBox();
    void loadBonds(string filename);
    void loadAngles(string filename);
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

        while(!in.eof())
        {
            bond = Bond();
            in >> bond.N >> bond.type >> bond.at1 >> bond.at2;
            if(bond.N >= 0)
                bonds.push_back(bond);
        }
    }
    in.close();

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
