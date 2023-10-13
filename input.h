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

using namespace std;

/**
 * @brief The BeadParam class
 * Bead LJ parameters
 */
class LJParam {
public:
    LJParam(){}
    LJParam(int type, double epsilon, double sigma, double cutoff) : type(type), epsilon(epsilon), sigma(sigma), cutoff(cutoff) {}

    int type=-1;
    double epsilon=1.0;
    double sigma=1.0;
    double cutoff=25.0;

    string toString()
    {
        stringstream ss;
        ss << type << " " << epsilon << " " << sigma << " " << cutoff << endl;
        return ss.str();
    }
};

class CosParam {
public:
    CosParam(){}
    CosParam(int type1, int type2, double epsilon, double start_dis, double range) : type1(type1), type2(type2), epsilon(epsilon), start_dis(start_dis), range(range){}

    int type1;
    int type2;
    double epsilon;
    double start_dis;
    double range;

    string toString()
    {
        stringstream ss;
        ss << type1 << " " << type2 << " " << epsilon << " " << start_dis << " " << range << endl;
        return ss.str();
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
 * @brief The Input class - Parameters for particle generation
 */
class Input{
public:
    Input() {}

    Output out; /// Output type - none, pdb, lammps_full, xyz

    int nano_type;

    int num_of_beads;
    myFloat scale = 0.0;
    int offset = 1;
    int seed=0;
    myFloat c;
    int num_lig;
    int chain_type=-1;
    int mol_tag=-1;
    int atom_type=1;
    bool center=false;
    int mtag_1=-1;
    int mtag_2=-1;
    Atom ivx = Atom(0.0, 0.0, 0.0);
    bool fit = false;
    Atom boxm = Atom(-1,-1,-1);
    Atom boxp = Atom(-1,-1,-1);
    Atom com_pos = Atom(0.0, 0.0, 0.0);
    Atom patch_1 = Atom(1,1,1,0);
    Atom patch_2 = Atom(1,1,1,0);
    string infile;

    vector<LJParam> bparam;
    vector<CosParam> cparam;
    vector<int> types;

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
            if( what.compare("Particle_type:") == 0 )   { ss >> nano_type; }
            if( what.compare("Output_type:") == 0 )     { ss >> out; }
            if( what.compare("Num_of_beads:") == 0 )    { ss >> num_of_beads; }
            if( what.compare("Scale:") == 0 )           { ss >> scale; }
            if( what.compare("Lammps_offset:") == 0 )   { ss >> offset; }
            if( what.compare("c:") == 0 )               { ss >> c; }
            if( what.compare("Number_of_ligands:") == 0 ) { ss >> num_lig; }
            if( what.compare("Chain_type:") == 0 ) { ss >> chain_type; }
            if( what.compare("Box:") == 0 )             { ss >> boxm.x >> boxp.x >> boxm.y >> boxp.y >> boxm.z >> boxp.z; }
            if( what.compare("Position_shift:") == 0 )  { ss >> com_pos.x >> com_pos.y >> com_pos.z; }
            if( what.compare("Load_file:") == 0 )  { ss >> infile; }
            if( what.compare("Center") == 0 )  { center=true; }
            if( what.compare("Fit") == 0 )  { fit=true; }
            if( what.compare("Align:") == 0 )  { ss >> mtag_1 >> mtag_2;  }
            if( what.compare("Impact_vector:") == 0 )  { ss >> ivx.x >> ivx.y >> ivx.z; }
            if( what.compare("Patch_1:") == 0 )  { ss >> patch_1.vx >> patch_1.x >> patch_1.vy >> patch_1.y >> patch_1.vz >> patch_1.z >> patch_1.type; }
            if( what.compare("Patch_2:") == 0 )  { ss >> patch_2.vx >> patch_2.x >> patch_2.vy >> patch_2.y >> patch_2.vz >> patch_2.z >> patch_2.type; }
            if( what.compare("Seed:") == 0 )  { ss >> seed; rng.seed(seed); }

            if( what.compare("Mol_tag:") == 0 ) { ss >> mol_tag; }
            if( what.compare("Atom_type:") == 0 ) { ss >> atom_type; }

            if( what.compare("Beads_lj/cut:") == 0 )
            {
                ss >> len;
                for(int i=0; i < len; ++i)
                {
                    ss.flush();
                    ss.clear();
                    getline(fs, line);
                    ss.str(line);
                    bparam.push_back(LJParam());
                    ss >> bparam.back().type >> bparam.back().epsilon >> bparam.back().sigma >> bparam.back().cutoff;
                }
                len=0;
            }

            if( what.compare("Cosatt:") == 0 )
            {
                ss >> len;
                for(int i=0; i < len; ++i)
                {
                    ss.flush();
                    ss.clear();
                    getline(fs, line);
                    ss.str(line);
                    cparam.push_back(CosParam());
                    ss >> cparam.back().type1 >> cparam.back().type2 >> cparam.back().epsilon >> cparam.back().start_dis >> cparam.back().range;
                }
            }

            what.clear();
        }
        fs.close();

        for(auto i : bparam)
        {
            types.push_back(i.type);
        }
        return true;
    }

    string toString()
    {
        stringstream ss;

        ss << "Particle_type: " << nano_type << endl;
        ss << "Output_type: " << out << endl;
        ss << "Number of beads: " << num_of_beads << endl;
        ss << "Scale: " << scale << endl;
        ss << "Offset: " << offset << endl;
        ss << "c: " << c << endl;
        ss << "Number of ligands: " << num_lig << endl;
        ss << "Box: (" << boxm.x << ", " << boxp.x << ", " << boxm.y << ", " << boxp.y << ", " << boxm.z << ", " << boxp.z << ")" << endl;
        ss << "Position: (" << com_pos.x << ", " << com_pos.y << ", " << com_pos.z << ")" << endl;
        ss << "Patch_1: (" << patch_1.x << "-" << patch_1.vx << ", " << patch_1.y << "-" << patch_1.vy << ", " << patch_1.z << "-" << patch_1.vz << ", " << patch_1.type << ")" << endl;
        ss << "Patch_1: (" << patch_2.x << ", " << patch_2.y << ", " << patch_2.z << ", " << patch_2.type << ")" << endl;

        ss << "Beads_lj/cut:" << endl;
        for(auto i : bparam)
        {
            ss << i.toString();
        }

        ss << "Cosatt:" << endl;
        for(auto i : cparam)
        {
            ss << i.toString();
        }

        return ss.str();
    }

    void clear()
    {
        nano_type=-1;
        out.clear();
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

        boxm=Atom(-1,-1,-1);
        boxp=Atom(-1,-1,-1);
        com_pos=Atom(0.0, 0.0, 0.0);
        patch_1=Atom(1,1,1,0);
        patch_2=Atom(1,1,1,0);
        ivx=Atom(0.0, 0.0, 0.0);

        infile.clear();
        bparam.clear();
        cparam.clear();
        types.clear();
    }

    bool is_mtag_12()
    {
        if( mtag_1 != -1 && mtag_2 != -1 )
            return false;
        return true;
    }

    bool is_mol_tag()
    {
        if(mol_tag == -1)
            return false;
        return true;
    }

    bool isCOM_pos()
    {
        if(com_pos.x == 0.0 && com_pos.y == 0.0 && com_pos.z == 0.0)
            return false;
        return true;
    }

    bool isScale()
    {
        if(scale == 0.0)
            return false;
        return true;
    }
};


#endif // INPUT_H
