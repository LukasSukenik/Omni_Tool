#ifndef PHAGE_MEM_CLOSURE_H
#define PHAGE_MEM_CLOSURE_H

#include "system_base.h"
#include "atom.h"

class Phage_Mem_Closure : public System_Base
{
public:
    inline static const string keyword = "Phage_Mem_Closure";
    const string name = "Phage_Mem_Closure";

    Phage_Mem_Closure() : System_Base("Phage_Mem_Closure") {}

    string help()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Phage_Mem_Closure" << endl;
        ss << "System_execute: Smart_convert_pdb" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate(Data& data)
    {
        data.in.param.validate_keyword("Load_file", "data.start");
    }

    void execute(Data& data)
    {
        validate(data);
        if(data.in.param["System_execute"].compare("Convert_pdb") == 0) { smart_convert_pdb(data); }
    }

private:

    /**
     * @brief smart_convert_pdb - Convert molecules atom.type, so that they can be easily distinguished using type only
     * - as input we need atom.type, atom.mol_tag, and bond information to disambiguate
     * @param data
     *
     * Molecules:
     * atom.type 1-2-2-3,   atom.mol_tag = 1 -> Normal Lipid  -> 1-2-2-3, mol_tag=1
     * atom.type 4-2-2-3,   atom.mol_tag = 1 -> Lipid II carrier -> 4-5-5-6, mol_tag=2
     * atom.type 5 6 7 8 9, atom.mol_tag = 2 -> 7 8 9 10 11, mol_tag=3
     * 5 - Rigid bacteriophage tail
     * 6 - Netted-slab - Body
     * 7 - Netted-slab - Interacting edges
     * 8 - Netted-slab - Lipid Anchor Head
     * 9 - Netted-slab - Lipid Anchor Tail
     *
     * Rest of molecules have incrementally larger mol_tag, because they are GCMC inserts
     *
     * 11-10-12 - BB_Free                     -> 13-12-14, mol_tag=4
     * 4-(b8)-13-10-12 - BB_Mem               -> 4-(b8)-16-15-17, mol_tag=5
     * 4-(b8)-13-10-14-(b10)-6 - BB_Mem_Slab  -> 4-(b8)-19-18-20-(b10)-8, mol_tag=6
     * 4-(b8)-13-10-14-(b12)-16 - BB_Mem_Slab -> 4-(b8)-19-18-20-(b12)-24
     * 15-10-14-(b10)-6 - BB_Slab             -> 22-21-23-(b10)-8, mol_tag=7
     * 16-10-14-(b10)-6 - BB_Slab             -> 24-21-23-(b10)-8
     * 15-10-14-(b12)-16 - BB_Slab            -> 22-21-23-(b12)-24
     * 16-10-14-(b12)-16 - BB_Slab            -> 24-21-23-(b12)-24
     *
     */
    void smart_convert_pdb(Data& data)
    {
        cerr << "Phage_Mem_Closure::smart_convert_pdb" << endl;

        Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        //Bonds& bonds = data.coll_bonds[   data.id_map[ data.in.p_int["ID"] ]   ];


        // loop over topo and change types without molecule check
        for(size_t i=0; i<topo.size(); ++i)
        {
            convert_BB_Slab(topo,i);           // 22-21-23 and 24-21-23
            convert_BB_Mem_Slab(topo,i);       // 19-18-20
            convert_BB_Mem(topo,i);            // 16-15-17
            convert_BB_Free(topo, i);          // 13-12-14
            convert_wall_phage_tail(topo, i);  // 7 8 9 10 11
        }
        for(size_t i=0; i<topo.size(); ++i)
        {
            convert_lipid_II_carrier(topo, i); // 4-5-5-6
        }
    }

    void convert_lipid_II_carrier(Atoms& topo, size_t i)
    {
        // atom.type 4-2-2-3,   atom.mol_tag = 1 -> Lipid II carrier -> 4-5-5-6, mol_tag=2
        if(topo[i].type == 4)
        {
            if(topo[i+1].type == 2 && topo[i+2].type == 2 && topo[i+3].type == 3)
            {
                topo[i+0].mol_tag=2;
                topo[i+1].mol_tag=2; topo[i+1].type = 5;
                topo[i+2].mol_tag=2; topo[i+2].type = 5;
                topo[i+3].mol_tag=2; topo[i+3].type = 6;
            }
            else
            {
                cerr << "Lipid II Carrier type change failed" << endl;
                exit(-1);
            }
        }
    }

    void convert_wall_phage_tail(Atoms& topo, size_t i)
    {
        // atom.type 5 6 7 8 9, atom.mol_tag = 2 -> 7 8 9 10 11, mol_tag=3
        if(topo[i].type == 9 && topo[i].mol_tag == 2) { topo[i].mol_tag=3; topo[i].type = 11; }
        if(topo[i].type == 8 && topo[i].mol_tag == 2) { topo[i].mol_tag=3; topo[i].type = 10; }
        if(topo[i].type == 7 && topo[i].mol_tag == 2) { topo[i].mol_tag=3; topo[i].type = 9; }
        if(topo[i].type == 6 && topo[i].mol_tag == 2) { topo[i].mol_tag=3; topo[i].type = 8; }
        if(topo[i].type == 5 && topo[i].mol_tag == 2) { topo[i].mol_tag=3; topo[i].type = 7; }
    }

    void convert_BB_Free(Atoms& topo, size_t i)
    {
        // BB_Free 11-10-12 -> 13-12-14, mol_tag=4
        if(topo[i].type == 11 &&
           topo[i+1].type == 10 &&
           topo[i+2].type == 12 &&
           (topo[i].mol_tag == topo[i+1].mol_tag && topo[i].mol_tag == topo[i+2].mol_tag) )
        {
            topo[i+0].mol_tag=4; topo[i+0].type=13; // 11
            topo[i+1].mol_tag=4; topo[i+1].type=12; // 10
            topo[i+2].mol_tag=4; topo[i+2].type=14; // 12
        }
    }

    void convert_BB_Mem(Atoms& topo, size_t i)
    {
        // BB_Mem 13-10-12 -> 16-15-17, mol_tag=5
        if(topo[i].type == 13 &&
           topo[i+1].type == 10 &&
           topo[i+2].type == 12 &&
           (topo[i].mol_tag == topo[i+1].mol_tag && topo[i].mol_tag == topo[i+2].mol_tag) )
        {
            topo[i+0].mol_tag=5; topo[i+0].type=16;
            topo[i+1].mol_tag=5; topo[i+1].type=15;
            topo[i+2].mol_tag=5; topo[i+2].type=17;
        }
    }

    void convert_BB_Mem_Slab(Atoms& topo, size_t i)
    {
        // BB_Mem_Slab 13-10-14 -> 19-18-20, mol_tag=6
        if(topo[i].type == 13 &&
           topo[i+1].type == 10 &&
           topo[i+2].type == 14 &&
           (topo[i].mol_tag == topo[i+1].mol_tag && topo[i].mol_tag == topo[i+2].mol_tag) )
        {
            topo[i+0].mol_tag=6; topo[i+0].type=19;
            topo[i+1].mol_tag=6; topo[i+1].type=18;
            topo[i+2].mol_tag=6; topo[i+2].type=20;
        }
    }

    void convert_BB_Slab(Atoms& topo, size_t i)
    {
        // BB_Slab 16-10-14 -> 24-21-23, mol_tag=7
        if(topo[i].type == 16 &&
           topo[i+1].type == 10 &&
           topo[i+2].type == 14 &&
           (topo[i].mol_tag == topo[i+1].mol_tag && topo[i].mol_tag == topo[i+2].mol_tag) )
        {
            topo[i+0].mol_tag=7; topo[i+0].type=24;
            topo[i+1].mol_tag=7; topo[i+1].type=21;
            topo[i+2].mol_tag=7; topo[i+2].type=23;
        }

        // BB_Slab 15-10-14 -> 22-21-23, mol_tag=7
        if(topo[i].type == 15 &&
           topo[i+1].type == 10 &&
           topo[i+2].type == 14 &&
           (topo[i].mol_tag == topo[i+1].mol_tag && topo[i].mol_tag == topo[i+2].mol_tag) )
        {
            topo[i+0].mol_tag=7; topo[i+0].type=22;
            topo[i+1].mol_tag=7; topo[i+1].type=21;
            topo[i+2].mol_tag=7; topo[i+2].type=23;
        }
    }

    void print_atom(Atoms& topo, size_t i)
    {
        cerr << topo[i].N << " " << topo[i].mol_tag << " " << topo[i].type << endl;
    }
};

#endif // PHAGE_MEM_CLOSURE_H
