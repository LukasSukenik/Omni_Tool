#ifndef IO_PDB_H
#define IO_PDB_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>

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

        string line;
        Atom part;
        std::fstream in;
        std::stringstream ss;

        in.open(in_file, std::fstream::in);

        if(in.is_open())
        {
            while(in.good())
            {
                getline(in, line);

                ss.str( line );

                // string indexed from 0
                if( string("ATOM  ").compare( line.substr(0, 6) ) == 0 ||  string("HETATM ").compare( line.substr(0, 6) ) == 0 )
                {
                    part.atom_serial_N = line.substr(6, 5);              // char  : 7-11: atom serial number (not decadic, contains letters)
                    part.atom_name = line.substr(12, 4);                 // char  : 13-16: atom name
                    part.res_name = line.substr(17, 3);                  // char  : 18-20: Residue name
                    part.chain_id = line[21];                            // char  : 22: chain identifier
                    part.res_seq_N = std::stoi(line.substr(22, 4));      // int   : 23-26: Residue sequence number,
                    part.code = line[26];                                // char  : 27:code for insertion of residues
                    part.pos.x = std::stof(line.substr(30, 8));          // float : 31-38
                    part.pos.y = std::stof(line.substr(38, 8));          // float : 39-46
                    part.pos.z = std::stof(line.substr(46, 8));          // float : 47-54
                    part.occupancy = std::stof(line.substr(54, 6));      // float : 55-60: Occupancy
                    part.temp_factor = std::stof(line.substr(60, 6));    // float : 61-66: Temperature factor
                    part.seg_id = line.substr(72, 4);                    // char : 73-76: Segment identifier (optional)
                    part.element = line.substr(76, 2);                   // char : 77-78 Element symbol
                    part.charge = line.substr(78, 2);                    // char : 79-80 Charge (optional)

                    beads.push_back(part);
                }

                ss.flush();
                ss.clear();


            }
            cerr << "  Added " << beads.size() << " beads" << endl;
        } else {
            cerr << "File " << in_file << " file not opened" << endl;
        }

        in.close();
    }

    void print_atom(Atom& atom)
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

    void print()
    {
        for (Atom& atom : beads)
        {
            print_atom(atom);
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

    Atoms beads;

    void print()
    {
        cout << beads.size() << "\nparticle\n";
        for (const Atom& atom : beads)
        {
            cout << "C" << atom.type <<  " " << atom.pos << "\n";
        }
    }
};

#endif // IO_PDB_H
