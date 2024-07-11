#ifndef IO_PDB_H
#define IO_PDB_H

#include <iostream>
#include <iomanip>

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

    void print(Atoms& all_beads)
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

        for (Atom& atom : all_beads) {
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
