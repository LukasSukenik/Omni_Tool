#ifndef IO_PDB_H
#define IO_PDB_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string.h>

#include "atom.h"

using namespace std;

/**
 * @brief The PDB class
 * Class for IO of files in PDB data format
 */
class IO_PDB
{
public:
    IO_PDB(){}

    Atoms beads;

    void load(string in_file)
    {
        if(in_file.empty())
        {
            return;
        }

        char str[256];
        Atom part;
        std::fstream in;
        std::stringstream ss;

        in.open(in_file, std::fstream::in);

        if(in.is_open())
        {
            // Load until "Atoms" occurs
            while(in.good() && strstr(str, "Atoms") == NULL)
            {
                in.getline(str, 256);
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
            cerr << "  Added " << beads.size() << " beads" << endl;
        } else {
            cerr << "File " << in_file << " file not opened" << endl;
        }

        in.close();
    }

    void print(Atoms& all_beads)
    {
        for (Atom& atom : all_beads)
        {
            cout << "ATOM"                                             // 1-4
                 << "  "                                               // 5-6: empty
                 << std::setw(5) << std::right << atom.atom_serial_N   // 7-11: atom serial number,
                 << " "                                                // 12: empty
                 << std::setw(4) << std::left  << atom.atom_name       // 13-16: atom name
                 << " "                                                // 17: empty
                 << std::setw(3) << std::right << atom.res_name        // 18-20: Residue name
                 << " "                                                // 21: empty
                 << std::setw(1)               << atom.chain_id        // 22: chain identifier
                 << std::setw(4) << std::right << atom.res_seq_N       // 23-26: Residue sequence number
                 << std::setw(1)               << atom.code            // 27:code for insertion of residues
                 << "   "                                              // 28-30:empty
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.pos.x // 31-38
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.pos.y // 39-46
                 << std::setw(8) << std::fixed << std::setprecision(3) << atom.pos.z // 47-54
                 << std::setw(6) << std::right << atom.occupancy       // 55-60: Occupancy
                 << std::setw(6) << std::right << atom.temp_factor     // 61-66: Temperature factor
                 << "      "                                           // 67-72: Empty
                 << std::setw(4) << std::left  << atom.seg_id          // 73-76: Segment identifier (optional)
                 << std::setw(2) << std::right << atom.element         // 77-78 Element symbol
                 << std::setw(2)               << atom.charge << "\n"; // 79-80 Charge (optional)
        }
    }

    void print_lammps_data(Atoms& all_beads)
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

        for (Atom& atom : all_beads)
        {
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
};

class IO_XYZ
{
public:
    IO_XYZ() {}

    void print(Atoms& all_beads)
    {
        cout << all_beads.size() << "\nparticle\n";
        for (const Atom& atom : all_beads)
        {
            cout << "C" << atom.type <<  " " << atom.pos << "\n";
        }
    }
};

#endif // IO_PDB_H
