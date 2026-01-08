#ifndef LIPID_H
#define LIPID_H

#include <vector>

#include "particle.h"
#include "atom.h"
#include "bond.h"

#define RECEPTOR_BOND_TYPE 4
#define RECEPTOR_ANGLE_TYPE 1

using namespace std;


/**
 * @brief The Lipid class - Deserno 4 bead lipid
 *
 * TODO:
 * - only works when generating into an empty system, implement particle and bond offseting pending
 * - move particle type setting based on leaflet into constructor
 * - set_input - differentiate WHZ and NO WHZ in input
 */
class Lipid : public Particle {
public:
    inline static const string keyword = "lipid";

    enum class Leaflet {upper, lower};

    ///
    /// default bead types for 4 bead model without WHZ - 1-2-2-3 || 3-2-2-1
    /// Default bead types for 4 bead model with WHZ - 1-2-2-3 || 6-5-5-4
    ///
    int type_head_upper_leaf = 1;
    int type_tail_upper_leaf = 2;
    int type_tail_2_upper_leaf = 3;
    int type_head_lower_leaf = 4;
    int type_tail_lower_leaf = 5;
    int type_tail_2_lower_leaf = 6;
    int type_receptor = 7;

    ///
    /// Default bond types
    ///
    int bond_type_tail_head = 1;
    int bond_type_tail_tail = 2;
    int bond_type_harmonic = 3;
    int bond_type_tail_end = 4;
    int bond_type_harmonic2 = 5;

    Atoms part;  // consecutive:  0:HEAD, 1:TAIL, 2:TAIL (the further one)
    Bonds bond;  // 0: between 1-2 TAIL_HEAD_BOND, 1: 2-3 TAIL_TAIL_BOND, 2: 1-3 HARMONIC_BOND, 3: 3-4 TAIL_END_BOND, 4: 2-4 HARMONIC_BOND2


    Lipid() : Particle("lipid") {}
    Lipid(Atom p1, Atom p2, Atom p3, Atom p4, int mol_N=0, Leaflet leaf=Leaflet::upper, bool is_receptor=false) : Particle("lipid")
    {
        // assumes atom.pos, atom.N, atom.mol_tag is set, atom.type is now set after in vesicle.h, TODO move inside constructor
        store_lipid_part(p1, p2, p3, p4);
        set_lipid_part_N(mol_N);
        set_bead_type(leaf, is_receptor);

        // set bond.N, bond.N, and bond.type
        bond.resize(5);
        set_bond_types();
        set_lipid_bond_N(mol_N);
    }
    Lipid(Tensor_xyz pos, Tensor_xyz dir, int mol_N=0, int mol_tag=1, Leaflet leaf=Leaflet::upper, bool is_receptor=false) : Particle("lipid")
    {
        gen_lipid(pos, dir, mol_N, mol_tag); set_bead_type(leaf, is_receptor);
    }

    /**
     * @brief generate
     * @param data
     * Deserno lipid in Lammps is defined by
     * - Atoms
     * -- atom.N
     * -- atom.pos -> set directly, or better by position of head and direction
     * -- atom.type -> type_head_upper_leaf, ...
     * -- atom.mol_tag -> set for convenience, not necessary but nice to have
     * - Bonds
     * -- bond.N,
     * -- bond.type,
     * -- bond.at1 and at2 ->
     */
    void generate( Data& data )
    {
        set_input(data);
        gen_lipid(Tensor_xyz(0,0,0), Tensor_xyz(0,0,1), 0, data.in.mol_tag);

        beads.insert(beads.end(), part.begin(), part.end());
        bonds.insert(bonds.end(), bond.begin(), bond.end());
    }

private: // Lipid Generation
    void set_input( Data& data )
    {
        type_head_upper_leaf = data.in.atom_type[0];
        type_tail_upper_leaf = data.in.atom_type[1];
        type_tail_2_upper_leaf = data.in.atom_type[2];
        type_head_lower_leaf = data.in.atom_type[3];
        type_tail_lower_leaf = data.in.atom_type[4];
        type_tail_2_lower_leaf = data.in.atom_type[5];
        type_receptor = data.in.atom_type[6];
    }

    void set_bond_types()
    {
        bond[0].type = bond_type_tail_head;
        bond[1].type = bond_type_tail_tail;
        bond[2].type = bond_type_harmonic;
        bond[3].type = bond_type_tail_end;
        bond[4].type = bond_type_harmonic2;
    }

    void set_bead_type(Leaflet leaf, bool is_receptor=false)
    {
        if(leaf == Leaflet::upper)
        {
            part[0].type = type_head_upper_leaf;
            part[1].type = type_tail_upper_leaf;
            part[2].type = type_tail_upper_leaf;
            part[3].type = type_tail_2_upper_leaf;
        }
        if(leaf == Leaflet::lower)
        {
            part[0].type = type_head_lower_leaf;
            part[1].type = type_tail_lower_leaf;
            part[2].type = type_tail_lower_leaf;
            part[3].type = type_tail_2_lower_leaf;
        }
        if(is_receptor)
        {
            part[0].type = type_receptor;
        }
    }



private: // Lipid Generation
    void gen_lipid(Tensor_xyz pos, Tensor_xyz dir, int mol_N=0, int mol_tag=1)
    {
        // Atom.N, Atom.pos, Atom.type, Atom.mol_tag
        Atom p1 = Atom(4*mol_N +1, pos, type_head_upper_leaf, mol_tag);
        Atom p2 = Atom(4*mol_N +2, pos + dir, type_tail_upper_leaf, mol_tag);
        Atom p3 = Atom(4*mol_N +3, pos + dir*2.0, type_tail_upper_leaf, mol_tag);
        Atom p4 = Atom(4*mol_N +4, pos + dir*3.0, type_tail_2_upper_leaf, mol_tag);

        // internal store withing the lipid class
        store_lipid_part(p1, p2, p3, p4);

        bond.resize(5);
        set_bond_types(); // bond.type
        set_lipid_bond_N(4*mol_N +1, 5*mol_N +1); // bond.N, bond.at1, bond.at2
    }

    void set_bonds();

    void store_lipid_part(Atom p1, Atom p2, Atom p3, Atom p4)
    {
        part.resize(4);

        part[0] = p1;
        part[1] = p2;
        part[2] = p3;
        part[3] = p4;
    }

    void set_lipid_part_N(int mol_N=0)
    {
        part[0].N = 4*mol_N +1;
        part[1].N = 4*mol_N +2;
        part[2].N = 4*mol_N +3;
        part[3].N = 4*mol_N +4;
    }

    void set_lipid_bond_N(int mol_N=0)
    {
        set_lipid_bond_N(4*mol_N +1, 5*mol_N +1);
    }

    void set_lipid_bond_N(int particle_N, int bond_N)
    {
        bond[0].N = bond_N;
        bond[1].N = bond_N+1;
        bond[2].N = bond_N+2;
        bond[3].N = bond_N+3;
        bond[4].N = bond_N+4;

        bond[0].at1 = particle_N;
        bond[0].at2 = particle_N+1;

        bond[1].at1 = particle_N+1;
        bond[1].at2 = particle_N+2;

        bond[2].at1 = particle_N;
        bond[2].at2 = particle_N+2;

        bond[3].at1 = particle_N+2;
        bond[3].at2 = particle_N+3;

        bond[4].at1 = particle_N+1;
        bond[4].at2 = particle_N+3;
    }


public:
    void report(Data& data)
    {
        if(beads.empty())
        {
            cerr << "No particle created, can't populate system" << endl;
            exit(-1);
        }
        else
        {
            cerr << "Population particle of " << part.size() << " atoms " << data.in.population.count << " times, coll_beads.size" << data.coll_beads.size() << endl;
            beads.clear(); // clear beads and bonds to generate new particles with random distribution
            bonds.clear(); // in lipid class the necessary particles and bodns already stored in part and bond
        }
    }

    int get_lipid_count(Data& data)
    {
        if(!data.coll_beads.empty())
        {
            int existing_lipid_count = data.coll_beads[0].size() / part.size();
            cerr << "lipid::populate coll_beads[0].size " << data.coll_beads[0].size() << "\n";
            return existing_lipid_count;
        }
        return 0;
    }

    void eval_attempts(int tries, int mol_count)
    {
        if(tries != 0 && tries % 10000 == 0)
            cerr << "Molecule " << mol_count << " trial " << tries << "\n";

        if(tries > 10*1000*1000) // with cell list 10M is feasible
        {
            cerr << "Particles::populate -> Can't generate any more particles, simulation box is full, tries" << tries << endl;
            exit(1);
        }
    }

    bool is_coll_overlap(vector<Atoms>& coll_Atoms)
    {
        if(!coll_Atoms.empty())
        {
            return coll_Atoms[0].is_overlap(part, 0.9, coll_cell_list.neighbors);
        }
        return false;
    }

    void populate(Data& data)
    {
        if( data.in.population.random )
        {
            report(data);

            Tensor_xyz random_pos;
            Tensor_xyz random_dir;
            int mol_tag = data.in.mol_tag;
            int existing_lipid_count = get_lipid_count(data);

            cell_list.init(data);
            coll_cell_list.init(data);

            if(!data.coll_beads.empty())
            {
                coll_cell_list.add(data.coll_beads[0]);
            }

            //
            // Copy, move, rotate
            //
            for(int i=0; i < data.in.population.count; ++i)
            {
                cerr << "lipid::populate lipid_N " << i+existing_lipid_count << "\n";

                random_pos = data.in.sim_box.get_random_pos();
                random_dir.randomUnitSphere();
                gen_lipid(random_pos, random_dir, i+existing_lipid_count, mol_tag); // stores new lipid into part

                cell_list.set_neighbors( part.get_center_of_mass().pos );
                coll_cell_list.set_neighbors( part.get_center_of_mass().pos );

                int tries=0;
                while( !data.in.sim_box.is_in_box(part) || beads.is_overlap(part, 0.9, cell_list.neighbors) || is_coll_overlap( data.coll_beads ) ) // beads.is_overlap(part, 0.9)
                {
                    random_pos = data.in.sim_box.get_random_pos();
                    random_dir.randomUnitSphere();
                    gen_lipid(random_pos, random_dir, i+existing_lipid_count, mol_tag);

                    cell_list.set_neighbors( part.get_center_of_mass().pos );
                    coll_cell_list.set_neighbors( part.get_center_of_mass().pos );

                    eval_attempts(tries, i);
                    ++tries;
                }

                cell_list.add(part, beads.size());
                beads.insert( beads.end(), part.begin(), part.end() );
                bonds.insert( bonds.end(), bond.begin(), bond.end() );

                cerr << "Generating population, particle " << i << " of " << data.in.population.count << ", successful on try " << tries << endl;

                //if(i==10) exit(1);
            }

            if(data.coll_beads.size() > 2)
            {
                cerr << "Lipid::populate - multiple collection overlap not implemented, max collections is 2" << endl;
            }

            cell_list.delete_all();
            coll_cell_list.delete_all();
        }
        cerr << "populate " << data.in.population << endl;
    }



    Atom get_direction()
    {
        Atom dir = part[0] - part[3]; // head - tail
        dir.normalise();
        return dir;
    }




};





class Lipids : public vector<Lipid>
{
public:
    bool update_positions(vector<Tensor_xyz>& frame)
    {
        if(4*size() != frame.size())
        {
            return false;
        }

        for(unsigned int i=0; i<size(); ++i) // loop over lipids
        {
            Lipid& lip = this->at(i);
            lip.part[0].pos = frame[4*i +0];
            lip.part[1].pos = frame[4*i +1];
            lip.part[2].pos = frame[4*i +2];
            lip.part[3].pos = frame[4*i +3];
        }

        return true;
    }

    Atoms get_Atoms()
    {
        Atoms at;
        for(Lipid& lip : (*this))
        {
            at.insert(at.end(), lip.part.begin(), lip.part.end());
        }
        return at;
    }

    void convert_receptors(int receptor_count)
    {
        int lipid_count = size();
        int count=0;
        int random=0;
        while(count < receptor_count)
        {
            random = (int)(ran() * lipid_count); // select random lipid
            if(random < size())
            {
                if(this->at(random).part[0].type == this->at(random).type_head_upper_leaf)
                {
                    this->at(random).part[0].type = this->at(random).type_receptor;
                    count++;
                }
                if(this->at(random).part[0].type == this->at(random).type_head_lower_leaf)
                {
                    this->at(random).part[0].type = this->at(random).type_receptor;
                    count++;
                }
            }
        }
    }
};




bool isTop(int i, vector<Atom>& particles, int head_upper_leaf=1, int head_lower_leaf = 4) {
    //
    // search particles, if another particle of type HEAD found under i, then i is top layer
    //
    cerr << "lipid.h :: isTop function broken, head_upper_leaf and head_lower_leaf changed to random default" << endl;

    for(unsigned int j=0; j<particles.size(); ++j) {
        if(particles[j].type == head_upper_leaf ||  particles[j].type == head_lower_leaf) {
            if( ( (particles[j].pos.x - particles[i].pos.x)*(particles[j].pos.x - particles[i].pos.x)
                  + (particles[j].pos.y - particles[i].pos.y)*(particles[j].pos.y - particles[i].pos.y)
                  + (particles[j].pos.z - particles[i].pos.z)*(particles[j].pos.z - particles[i].pos.z) ) < 36.0 ) {
                if( (particles[j].pos.z - 4.0 > particles[i].pos.z) ) {
                    return true;
                }
            }
        }
    }
    return false;
}

#endif // LIPID_H
