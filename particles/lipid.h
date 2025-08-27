#ifndef LIPID_H
#define LIPID_H

#include <vector>

#include "particle.h"
#include "atom.h"
#include "bond.h"



#define RECEPTOR_BOND_TYPE 4
#define RECEPTOR_ANGLE_TYPE 1

#define TAIL_HEAD_BOND 1
#define TAIL_TAIL_BOND 2
#define HARMONIC_BOND 3
#define TAIL_END_BOND 4
#define HARMONIC_BOND2 5

using namespace std;

class Tensor_xyz_integer
{
public:
    Tensor_xyz_integer(){}
    int x,y,z;
};

class Cell_List : public vector<vector<vector< vector<int >* >>> // 3D spatial indices + list of index to coll_beads
{
public:
    Cell_List(){}

    Tensor_xyz_integer number_of_cells;
    Tensor_xyz cell_size; // for lipids at leat 4.0

    double xlo = 0.0; // for spatial indices
    double xhi = 0.0;
    double ylo = 0.0;
    double yhi = 0.0;
    double zlo = 0.0;
    double zhi = 0.0;

    vector< vector<int>* > neighbors;
    vector<int>* empty_list;

    void show()
    {
        cerr << "xlo:" << xlo << " xhi:" << xhi << endl;
        cerr << "ylo:" << ylo << " yhi:" << yhi << endl;
        cerr << "zlo:" << zlo << " zhi:" << zhi << endl;
        cerr << "number of cells: " << number_of_cells.x << ", " << number_of_cells.y << ", " << number_of_cells.z << endl;
        cerr << "cell_size: " << cell_size.x << ", " << cell_size.y << ", " << cell_size.x << endl;
    }


    void init_cell_list(Data& data)
    {
        xlo = data.in.sim_box.xlo;
        xhi = data.in.sim_box.xhi;
        ylo = data.in.sim_box.ylo;
        yhi = data.in.sim_box.yhi;
        zlo = data.in.sim_box.zlo;
        zhi = data.in.sim_box.zhi;

        number_of_cells.x = floor( (xhi - xlo)/4.0 );
        number_of_cells.y = floor( (yhi - ylo)/4.0 );
        number_of_cells.z = floor( (zhi - zlo)/4.0 );

        cell_size.x = (xhi - xlo) / number_of_cells.x;
        cell_size.y = (xhi - xlo) / number_of_cells.x;
        cell_size.z = (xhi - xlo) / number_of_cells.x;

        empty_list = new vector<int>();

        neighbors.resize(27);
        neighbors.assign(27, empty_list);

        (*this).resize( number_of_cells.x );

        for(int i=0; i<number_of_cells.x; ++i)
        {
            (*this)[i].resize(number_of_cells.y);

            for(int j=0; j<number_of_cells.y; ++j)
            {
                (*this)[i][j].resize(number_of_cells.z);
                for(int k=0; k<number_of_cells.z; ++k)
                {
                    (*this)[i][j][k] = new vector<int >();
                    (*this)[i][j][k]->reserve(100);
                }
            }
        }
        cerr << "Cell_List::init_cell_list" << endl;
        show();
    }

    void delete_cell_list()
    {
        int number_of_cells = this->size();
        for(int i=0; i<number_of_cells; ++i)
        {
            for(int j=0; j<number_of_cells; ++j)
            {
                for(int k=0; k<number_of_cells; ++k)
                {
                    delete (*this)[i][j][k];
                }
            }
        }
        delete empty_list;
    }

    int get_spatial_X(Tensor_xyz pos)
    {
        return floor(floor(pos.x - xlo) / cell_size.x);
    }

    int get_spatial_Y(Tensor_xyz pos)
    {
        return floor(floor(pos.y - ylo) / cell_size.y);
    }

    int get_spatial_Z(Tensor_xyz pos)
    {
        return floor(floor(pos.z - zlo) / cell_size.z);
    }

    void add_lipid_to_cell_list(Atoms& part, int offset)
    {
        for(unsigned int i=0; i<part.size(); ++i)
        {
            (*this)[get_spatial_X(part[i].pos)]
                     [get_spatial_Y(part[i].pos)]
                     [get_spatial_Z(part[i].pos)]->push_back(offset + i);
            //cerr << "Cell_List::add_lipid_cell_list [" << get_spatial_index(part[0].pos.x) << ", " << get_spatial_index(part[0].pos.y) << ", " << get_spatial_index(part[0].pos.z) << "] = " << beads.size() + i << endl;
        }
    }

    void set_neighbors(Tensor_xyz pos)
    {
        //cerr << "Cell_List::set_neighbors start" << endl;
        int x = get_spatial_X(pos);
        int y = get_spatial_Y(pos);
        int z = get_spatial_Z(pos);
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                for(int k=0; k<3; ++k)
                {
                    if(x+i-1 < 0 || y+j-1 < 0 || z+k-1 < 0 || x+i-1 >= number_of_cells.x || y+j-1 >= number_of_cells.y || z+k-1 >= number_of_cells.z )
                    {
                        neighbors[(i*9) +(j*3) + k] = empty_list;
                    }
                    else
                    {
                        neighbors[(i*9) +(j*3) + k] = (*this)[x+i-1][y+j-1][z+k-1];
                    }
                }
            }
        }
        //cerr << "Cell_List::set_neighbors end" << endl;
        //show_neighbors();
    }

    void show_neighbors()
    {
        int sum=0;
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                for(int k=0; k<3; ++k)
                {
                    cerr << "Neighbor:" << i << ", " << j << ", " << k << " : " << neighbors[(i*9) +(j*3) + k]->size() << endl;
                    sum += neighbors[(i*9) +(j*3) + k]->size();
                }
            }
        }
        cerr << "Neighbor sum: " << sum << endl;
    }

    bool validate_neighbors(Atom test, Atoms& beads)
    {
        double max_dist_SQ = 4.0 * ( cell_size.x*cell_size.x + cell_size.y*cell_size.y + cell_size.z*cell_size.z );

        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                for(int k=0; k<3; ++k)
                {
                    for(int index : (*neighbors[(i*9) +(j*3) + k]) )
                    {
                        if(test.distSQ(beads[index]) > max_dist_SQ)
                        {
                            cerr << "validate_neighbors::Maximal distance exceeded: " << "[" << i << ", " << j << ", " << k << "] : " << sqrt(test.distSQ(beads[index])) << " > " << sqrt(max_dist_SQ) << endl;
                            cerr << "    test[" << test.pos.x << ", " << test.pos.y << ", " << test.pos.z << "] vs beads[" << index << "]:[" << beads[index].pos.x << ", " << beads[index].pos.y << ", " << beads[index].pos.z << "]" << endl;
                            exit(1);
                        }
                    }
                }
            }
        }
        return true;
    }
};

class Lipid : public Particle {
public:
    inline static const string keyword = "lipid";

    enum class Leaflet {upper, lower};

    int type_head = 1;
    int type_tail_1 = 2;
    int type_tail_2 = 3;

    static const int head_upper_leaf = 1;
    static const int tail_upper_leaf = 2;
    static const int tail_2_upper_leaf = 3;
    static const int head_lower_leaf = 4;
    static const int tail_lower_leaf = 5;
    static const int tail_2_lower_leaf = 6;
    static const int receptor = 7;

    Atoms part;  // consecutive:  0:HEAD, 1:TAIL, 2:TAIL (the further one)
    Bonds bond;      // 0: between 1-2 TAIL_HEAD_BOND, 1: 2-3 TAIL_TAIL_BOND, 2: 1-3 HARMONIC_BOND, 3: 3-4 TAIL_END_BOND, 4: 2-4 HARMONIC_BOND2

    Cell_List cell_list;


    Lipid() : Particle("lipid") {}
    Lipid(Atom p1, Atom p2, Atom p3, Atom p4) : Particle("lipid")
    {
        set_lipid(p1, p2, p3, p4);
    }


    void generate( Data& data )
    {
        set_input(data);
        set_lipid(0, Tensor_xyz(0,0,0), Tensor_xyz(0,0,1), data.in.mol_tag);

        beads.insert(beads.end(), part.begin(), part.end());
        bonds.insert(bonds.end(), bond.begin(), bond.end());
    }

    void set_input( Data& data )
    {
        type_head = data.in.atom_type[0];
        type_tail_1 = data.in.atom_type[1];
        type_tail_2 = data.in.atom_type[2];
    }

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




    void populate(Data& data)
    {
        if( data.in.population.random )
        {
            report(data);

            Tensor_xyz random_pos;
            Tensor_xyz random_dir;
            int mol_tag = data.in.mol_tag;
            int existing_lipid_count = get_lipid_count(data);

            cell_list.init_cell_list(data);

            //
            // Copy, move, rotate
            //
            for(int i=0; i < data.in.population.count; ++i)
            {
                random_pos = data.in.sim_box.get_random_pos();
                random_dir.randomUnitSphere();
                set_lipid(i+existing_lipid_count, random_pos, random_dir, mol_tag); // stores new lipid into part
                cell_list.set_neighbors( part.get_center_of_mass().pos );

                cerr << "lipid::populate lipid_N " << i+existing_lipid_count << "\n";

                int tries=0;
                bool is_coll_overlap = false;
                if(!data.coll_beads.empty())
                {
                    is_coll_overlap = data.coll_beads[0].is_overlap(part, 0.9);
                }
                while( !data.in.sim_box.is_in_box(part) || beads.is_overlap(part, 0.9, cell_list.neighbors) || is_coll_overlap ) // beads.is_overlap(part, 0.9)
                {
                    random_pos = data.in.sim_box.get_random_pos();
                    random_dir.randomUnitSphere();
                    set_lipid(i+existing_lipid_count, random_pos, random_dir, mol_tag);
                    cell_list.set_neighbors( part.get_center_of_mass().pos );

                    if(tries != 0 && tries % 10000 == 0)
                        cerr << "Molecule " << i << " trial " << tries << "\n";

                    if(tries > 1000000)
                    {
                        cerr << "Particles::populate -> Can't generate any more particles, simulation box is full, tries" << tries << endl;
                        exit(1);
                    }
                    ++tries;

                    if(!data.coll_beads.empty())
                    {
                        is_coll_overlap = data.coll_beads[0].is_overlap(part, 0.9);
                    }
                }

                cell_list.add_lipid_to_cell_list(part, beads.size());
                beads.insert( beads.end(), part.begin(), part.end() );
                bonds.insert( bonds.end(), bond.begin(), bond.end() );

                cerr << "Generating population, particle " << i << " of " << data.in.population.count << ", successful on try " << tries << endl;

                //if(i==10) exit(1);
            }

            if(data.coll_beads.size() > 2)
            {
                cerr << "Lipid::populate - multiple collection overlap not implemented, max collections is 2" << endl;
            }

            cell_list.delete_cell_list();
        }
        cerr << "populate " << data.in.population << endl;
    }

    void set_lipid(int mol_N, Tensor_xyz pos, Tensor_xyz dir, int mol_tag=1)
    {
        Atom p1 = Atom(4*mol_N +1, pos, type_head, mol_tag);
        Atom p2 = Atom(4*mol_N +2, pos + dir, type_tail_1, mol_tag);
        Atom p3 = Atom(4*mol_N +3, pos + dir*2.0, type_tail_1, mol_tag);
        Atom p4 = Atom(4*mol_N +4, pos + dir*3.0, type_tail_2, mol_tag);

        set_lipid(p1, p2, p3, p4, 4*mol_N +1, 5*mol_N +1);
    }


    void set_lipid(Atom p1, Atom p2, Atom p3, Atom p4, int pN=1, int bN=1)
    {
        part.resize(4);
        bond.resize(5);

        part[0] = p1;
        part[1] = p2;
        part[2] = p3;
        part[3] = p4;

        changeN_part(pN);
        changeN_bond(pN, bN);

        bond[0].type = TAIL_HEAD_BOND;
        bond[1].type = TAIL_TAIL_BOND;
        bond[2].type = HARMONIC_BOND;
        bond[3].type = TAIL_END_BOND;
        bond[4].type = HARMONIC_BOND2;
    }

    Atom get_direction()
    {
        Atom dir = part[0] - part[3]; // head - tail
        dir.normalise();
        return dir;
    }

    void changeN_part(int N) {
        part[0].N = N;
        part[1].N = N+1;
        part[2].N = N+2;
        part[3].N = N+3;
    }

    void changeN_bond(int pN, int bN) {
        bond[0].N = bN;
        bond[1].N = bN+1;
        bond[2].N = bN+2;
        bond[3].N = bN+3;
        bond[4].N = bN+4;

        bond[0].at1 = pN;
        bond[0].at2 = pN+1;

        bond[1].at1 = pN+1;
        bond[1].at2 = pN+2;

        bond[2].at1 = pN;
        bond[2].at2 = pN+2;

        bond[3].at1 = pN+2;
        bond[3].at2 = pN+3;

        bond[4].at1 = pN+1;
        bond[4].at2 = pN+3;
    }

    void set_bead_type(Leaflet leaf)
    {
        if(leaf == Leaflet::upper)
        {
            part[0].type = head_upper_leaf;
            part[1].type = tail_upper_leaf;
            part[2].type = tail_upper_leaf;
            part[3].type = tail_2_upper_leaf;
        }
        if(leaf == Leaflet::lower)
        {
            part[0].type = head_lower_leaf;
            part[1].type = tail_lower_leaf;
            part[2].type = tail_lower_leaf;
            part[3].type = tail_2_lower_leaf;
        }
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
};




bool isTop(int i, vector<Atom>& particles) {
    //
    // search particles, if another particle of type HEAD found under i, then i is top layer
    //
    for(unsigned int j=0; j<particles.size(); ++j) {
        if(particles[j].type == Lipid::head_upper_leaf ||  particles[j].type == Lipid::head_lower_leaf) {
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
