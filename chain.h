#ifndef CHAIN_H
#define CHAIN_H

#include "particle.h"

bool sort_by_dist(Atom& i, Atom& j)
{
    return i.temp_dist < j.temp_dist;
}


/**
 * @brief The Chain class
 * INPUT FILE:
 * Num_of_beads
 */
class Chain : public Particle
{
public:
    Chain() : Particle("Dodecahedron chain") {}
    Chain(string name) : Particle(name) {}

    const bool angle_pot = false;
    int chain_size= 3705;
    int core_size = 0;
    int chain_N = -1;
    const int moltype = 12;
    const double bead_size = 12.0;

    enum chain_type: unsigned short { BASIC=0, SOLID_CORE=1, INTERACTIVE_N=2 };
    chain_type chain_type = BASIC;

    void generate( Data& data )
    {
        load_params(data);
        erase_all_data(data);
        int type = data.all_sigma_size;
        set_FF(data);

        // Bond parameters
        double bond_size = 11.0;
        int offset = data.all_beads.size()+1;
        int bondOffset = data.all_bonds.size()+1;
        int bond_type = 86;

        //
        // Starting point, First chain BEAD
        //
        Atom next = get_COM(data);
        next.N = data.all_beads.size()+1;
        next.type = type+1;
        next.mol_tag = moltype;
        beads.push_back(next);

        //
        // Generating Chain
        //
        sigma[type][type] = 6.5; // needs to be below 11
        bool clashExist = true;
        int tick=0;
        for(int i=1; i<chain_size; ++i) {
            clashExist=true;
            while( clashExist ) {
                if(tick % 10000 == 0)
                {
                    cerr << "Try:" << tick << ", chain bead:" << i << endl;
                }
                if(tick > 100000)
                {
                    cerr << "Error, chain not generated" << endl;
                    exit(0);
                }

                // Generate new chain bead
                next.randomUnitSphere();
                next = next*bond_size + beads.back(); // move atom to last generated atom (+ convert * random unit dist)

                clashExist = ( clash(next, beads) || clash(next, data.all_beads) );

                //
                // Add bead to structure
                //
                if( !clashExist )
                {
                   next.mol_tag = moltype;
                   if( angle_pot ) {
                       next.type=type+2;
                   } else {
                       next.type=type+1;
                   }
                   next.N = beads.back().N + 1;
                   beads.push_back(next);
                   bonds.push_back( Bond(bonds.size() + bondOffset, bond_type, i+offset, i+offset-1, bond_size ) );

                   if( angle_pot ) {
                       //                                       current      -1        -2
                       angles.push_back( Angle(angles.size(), 1, i+offset, i+offset-1, i+offset-2 ) );
                   }
                } else {
                    ++tick;
                }
            }
        }

        sigma[type][type] = bead_size;
        cerr << beads.size() << endl;

        if(chain_type == SOLID_CORE)
        {
            Atom com = get_COM(data);
            vector<Atom> temp;
            int count=0;

            for(Atom& i : beads)
            {
                i.temp_dist = i.dist(com);
                temp.push_back(i);
            }
            sort(temp.begin(), temp.end(), sort_by_dist);

            for(Atom& i : temp)
            {
                if(count < core_size)
                for(Atom& j : beads)
                {
                    if(i.N == j.N)
                        j.type = type+2;
                }
                count++;
            }
        }
        if(chain_type == INTERACTIVE_N)
        {
            for(Atom& i : beads)
            {
                if(i.N % chain_N == 0)
                {
                    i.type = type+2;
                }
            }
        }
    }

    void load_params( Data& data )
    {
        chain_size = data.in.num_of_beads;
        core_size = data.in.num_lig;
        chain_N = data.in.num_lig;
        if(data.in.chain_type == 0)
            chain_type = BASIC;
        if(data.in.chain_type == 1)
            chain_type = SOLID_CORE;
        if(data.in.chain_type == 2)
            chain_type = INTERACTIVE_N;
    }

    void erase_all_data( Data& data )
    {
        // remove chain beads and bonds - type 12 and 86
        data.all_beads.erase( data.all_beads.begin()+get_num_capsid_beads(data), data.all_beads.end() );
        data.all_bonds.erase( data.all_bonds.begin()+get_first_chain_bond(data), data.all_bonds.end() );
    }

    int get_num_capsid_beads(Data& data)
    {
        int count = 0;
        for(Atom& i : data.all_beads)
        {
            if(i.type != 12)
                count++;
        }
        return count;
    }

    int get_first_chain_bond(Data& data)
    {
        int count=0;
        for(Bond& i : data.all_bonds)
        {
            if(i.type == 86)
                return count;
            count++;
        }
        return -1;
    }

    void set_FF( Data& data )
    {
        // Force field stuff - set sigma
        int type = data.all_sigma_size;
        loadFF(data);
        ++sigma_size;
        sigma[type][type] = bead_size;

        if( chain_type == SOLID_CORE || chain_type == INTERACTIVE_N )
        {
            ++sigma_size;
            sigma[type+1][type+1] = bead_size;
            mixing_rules(type+1);
        }

        mixing_rules(type);
        sigma_cosatt[11][11] = true;

        cerr << "Chain type:" << type << endl;
        printSigma(type);
    }

    void loadFF( Data& data )
    {
        sigma_size = data.all_sigma_size;
        for(int i=0; i<data.all_sigma_size; ++i)
        {
            for(int j=0; j<data.all_sigma_size; ++j)
            {
                sigma[i][j] = data.all_sigma[i][j];
                sigma_cosatt[i][j] = data.all_sigma_cosatt[i][j];
            }
        }
    }

    Atom get_COM(Data& data)
    {
        int count=0;
        Atom cm;
        for(Atom& item : data.all_beads)
        {
            cm += item;
            ++count;
        }
        cm *= 1.0/count;
        return cm;
    }


    bool clash(Atom& a, vector<Atom>& others)
    {
        bool clashExist = false;
        for( auto& o : others)
        {
            //cerr << a.dist(o) << " " << sigma[a.type-1][o.type-1] << " " << beads.size() << endl;
            if( a.dist(o) < sigma[a.type-1][o.type-1] ) // type inde in memory from 0, but lammps need types to start at 1
            {
                clashExist = true;
            }
        }
        return clashExist;
    }

    void mixing_rules(int i) {
        for(int j = 0; j< sigma_size; ++j) {
            if(i != j) {
                sigma[i][j] = 0.5*(sigma[i][i] + sigma[j][j]);
                sigma[j][i] = 0.5*(sigma[i][i] + sigma[j][j]);
            }
        }

    }
};

#endif // CHAIN_H
