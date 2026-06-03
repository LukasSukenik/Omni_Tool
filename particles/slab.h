#ifndef SLAB_H
#define SLAB_H

#include "particle.h"

class Slab1D : public Particle
{
public:
    inline static const string keyword = "slab1D";

    Slab1D() : Particle() {}

    string help()
    {
        stringstream ss;
        ss << "Particle_type: slab1D\n";
        return ss.str();
    }

    void generate( Data& data )
    {
        validate_inputs(data);
        int row = sqrt(data.in.p_int["Number_of_beads"]);
        int bead_size = data.in.p_float["Scale"];
        double factor = 1.0; // length were we generate beads


        for(int i=0; i<row; ++i) {
            for(int j=0; j<row; ++j) {
                for(int k=0; k<1; ++k) {
                     beads.push_back(Atom( factor*i, factor*k, factor*j, 0));
                }
            }
        }

        sigma_size = 1;
        sigma[0][0] = bead_size*1.0;
    }

private:
    void validate_inputs( Data& data )
    {
        if( !data.in.p_float.contains("Scale") )
        {
            cerr << "Missing keyword; Scale: 1.0" << endl;
            exit(-1);
        }
    }
};




class Slab : public Particle
{
public:
    inline static const string keyword = "slab";
    const string name = "slab";

    Slab() : Particle("slab") {}

    string help()
    {
        stringstream ss;
        ss << "Particle_type: slab\n";
        return ss.str();
    }

    void generate( Data& data )
    {
        validate_inputs(data);
        double factor = data.in.p_float["Scale"];

        sigma[0][0] = 5*factor;
        sigma[0][1] = 3*factor;
        sigma[1][1] = factor;
        sigma_size = 2;
        factor /= 4.0;

        vector<Atom> bVec;
        int size = data.in.p_int["Number_of_beads"];
        for(int a=0; a<size; ++a) {
            for(int b=0; b<size; ++b) {
                for(int c=0; c<size; ++c) {
                    if((a-size/2.0)*(a-size/2.0) + (b-size/2.0)*(b-size/2.0) + (c-size/2.0)*(c-size/2.0) < size*size*0.25)
                        bVec.push_back(Atom(a*factor*2.0 - size/2.0*factor, b*factor*2.0 - size/2.0*factor, c*factor*2.0 - size/2.0*factor, 1, 1));
                }
            }
        }

        Atom move = Atom( 14, 14, 14, 1, 1);
        for(Atom& item : bVec)
            item += move;

        beads.insert(this->beads.end(), bVec.begin(), bVec.end());
    }

private:
    void validate_inputs( Data& data )
    {
        if( !data.in.p_float.contains("Scale") )
        {
            cerr << "Missing keyword; Scale: 1.0" << endl;
            exit(-1);
        }
    }
};




class NettedSlab : public Particle
{
public:
    inline static const string keyword = "netted_slab";
    const string name = "netted_slab";

    NettedSlab() : Particle("netted_slab") {}

    string help()
    {
        stringstream ss;

        ss << "Particle_type: netted_slab\n";
        ss << "Output_type: lammps_full\n";
        ss << "Position_shift: 0 0 0\n";
        ss << "Scale: 1.0\n";
        ss << "Number_of_beads: 10\n";
        ss << "Mol_tag: 1\n";
        ss << "Atom_type: 1 2\n";
        ss << "Atom_type: 1->body 2->edges\n";

        return ss.str();
    }

    void validate_inputs( Data& data )
    {
        data.in.p_tensor.validate_keyword("Position_shift", "0 0 0");
        data.in.p_float.validate_keyword("Scale", "1.0");
        data.in.p_int.validate_keyword("Number_of_beads", "100");
        data.in.p_int.validate_keyword("Mol_tag", "1");
        data.in.p_vec_int.validate_keyword("Atom_type", "1 2");
        data.in.p_vec_int.validate_keyword("Bond_type", "1 2");
    }

    void generate( Data& data )
    {
        validate_inputs(data);

        Atoms slab;
        int row = sqrt(data.in.p_int["Number_of_beads"]);
        int N = 1; // offset handled by make persistent
        bool generated = false;

        for(int i=0; i<row; ++i) {
            for(int j=0; j<row; ++j) {
                for(int k=0; k<1; ++k) {
                    if( j==0 || i == 0 || i == row-1 || j == row-1 ) // edges
                    {
                        if(data.in.p_vec_int["Atom_type"].size() == 2)
                        {
                            slab.push_back(Atom(N, Tensor_xyz(i, j, k), data.in.p_vec_int["Atom_type"][1], data.in.p_int["Mol_tag"]));
                            generated = true;
                        }
                        if(data.in.p_vec_int["Atom_type"].size() == 5)
                        {
                            if(i == 0)     slab.push_back(Atom(N, Tensor_xyz(i, j, k), data.in.p_vec_int["Atom_type"][1], data.in.p_int["Mol_tag"]));
                            if(i == row-1) slab.push_back(Atom(N, Tensor_xyz(i, j, k), data.in.p_vec_int["Atom_type"][2], data.in.p_int["Mol_tag"]));
                            if(j == 0     && i != 0 && i != row-1) slab.push_back(Atom(N, Tensor_xyz(i, j, k), data.in.p_vec_int["Atom_type"][3], data.in.p_int["Mol_tag"]));
                            if(j == row-1 && i != 0 && i != row-1) slab.push_back(Atom(N, Tensor_xyz(i, j, k), data.in.p_vec_int["Atom_type"][4], data.in.p_int["Mol_tag"]));
                            generated = true;
                        }
                    }
                    else // body
                        slab.push_back(Atom(N, Tensor_xyz(i, j, k), data.in.p_vec_int["Atom_type"][0], data.in.p_int["Mol_tag"]));
                    ++N;
                }
            }
        }

        if(!generated)
        {
            cerr << "Missing edge atom type, Specify edge type, aka Atom_type: 1 2" << endl;
            exit(-1);
        }

        double radius = data.in.p_float["Radius"];
        double x,y;
        for(Atom& a : slab)
        {
            x = a.pos.x + (1.0/data.in.p_float["Scale"]) * data.in.p_tensor["Position_shift"].x;
            y = a.pos.y + (1.0/data.in.p_float["Scale"]) * data.in.p_tensor["Position_shift"].y;
            if(x*x + y*y < radius*radius)
            {
                a.type = 9;
                a.mol_tag = 1;
            }
        }

        gen_bonds(row, slab, data.in.p_vec_int["Atom_type"], data.in.p_vec_int["Bond_type"]);

        bonds.offset(get_coll_bond_N(data), get_coll_part_N(data));
        beads.insert(beads.end(), slab.begin(), slab.end());
    }

private:

    void mixing_rules() {
        for(int i = 0; i< sigma_size; ++i) {
            for(int j = 0; j< sigma_size; ++j) {
                if(i != j)
                    sigma[i][j] = 0.5*(sigma[i][i] + sigma[j][j]);
            }
        }
    }

    void gen_bonds(int row, Atoms& slab, vector<int> type_gen, vector<int> bond_type)
    {
        int part_ID;
        for(int i=0; i<row; ++i)
        {
            for(int j=0; j<row; ++j)
            {
                part_ID = i*row + j;

                if(is_bond_type(slab[part_ID].type, slab[(i+1)*row + j  ].type, type_gen) && i < row-1) bonds.push_back(    Bond(bonds.size()+1, bond_type[0], slab[(i+1)*row + j  ].N, slab[part_ID].N)    );
                if(is_bond_type(slab[part_ID].type, slab[(i  )*row + j+1].type, type_gen) && j < row-1) bonds.push_back(    Bond(bonds.size()+1, bond_type[0], slab[(i  )*row + j+1].N, slab[part_ID].N)    );

                if(bond_type.size() >= 2)
                {
                    if(is_bond_type(slab[part_ID].type, slab[(i+2)*row + j  ].type, type_gen) && i < row-2) bonds.push_back(    Bond(bonds.size()+1, bond_type[1], slab[(i+2)*row + j  ].N, slab[part_ID].N)    );
                    if(is_bond_type(slab[part_ID].type, slab[(i  )*row + j+2].type, type_gen) && j < row-2) bonds.push_back(    Bond(bonds.size()+1, bond_type[1], slab[(i  )*row + j+2].N, slab[part_ID].N)    );
                }
            }
        }
    }

    bool is_bond_type(int type, int type_next, vector<int> type_gen)
    {
        for(int tg : type_gen)
        {
            if(tg == type || tg == type_next)
                return true;
        }
        return false;
    }
};

#endif // SLAB_H
