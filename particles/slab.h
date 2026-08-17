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
        ss << "Output_type: lammps_full\n";
        ss << "Position_shift: 0 0 0\n";
        ss << "Scale: 1.0\n";
        ss << "Slab_size: 10 10 2\n";
        ss << "Mol_tag: 1\n";
        ss << "Atom_type: 1\n";
        ss << "Atom_type: (3 options):\n";
        ss << "define single type for body: 1\n";
        ss << "define 2 types: 1->body 2->edges\n";
        ss << "define 5 types: 1->body 2, 3, 4, 5->for each edge\n";

        return ss.str();
    }

    void validate_inputs( Data& data )
    {
        data.in.p_tensor.validate_keyword("Position_shift", "0 0 0");
        data.in.p_float.validate_keyword("Scale", "1.0");
        data.in.p_int.validate_keyword("Mol_tag", "1");
        data.in.p_vec_int.validate_keyword("Atom_type", "1 2");
        data.in.p_vec_int.validate_keyword("Slab_size", "10 10 2");

        if(data.in.p_vec_int["Slab_size"].size() != 3)
        {
            cerr << "Slab_size need to define 3 numbers. Slab_size: 10 10 2" << endl;
            exit(-1);
        }
    }

    void generate( Data& data )
    {
        validate_inputs(data);

        Atoms slab;
        int x_size = data.in.p_vec_int["Slab_size"][0];
        int y_size = data.in.p_vec_int["Slab_size"][1];
        int z_size = data.in.p_vec_int["Slab_size"][2];

        int N = 1; // offset handled by make persistent
        Tensor_xyz pos;
        int m_tag = data.in.p_int["Mol_tag"];
        int a_type = -1;

        int rr = 1;
        vector<int> x = {0+rr, 0+rr,     0+rr,        x_size/2, x_size/2,    x_size-rr-1, x_size-rr-1, x_size-rr-1};
        vector<int> y = {0+rr, y_size/2, y_size-rr-1, 0+rr,     y_size-rr-1, 0+rr,        y_size/2,    y_size-rr-1};

        for(int k=0; k<z_size; ++k)
        {
            for(int i=0; i<x_size; ++i)
            {
                for(int j=0; j<y_size; ++j)
                {
                    pos = Tensor_xyz(i, j, k);
                    a_type = data.in.p_vec_int["Atom_type"][0];

                    if(k==0) // ground floor
                    {
                        if(data.in.p_vec_int["Atom_type"].size() == 2) // 2 atom type have to be specified edges type
                        {
                            if(data.in.param["Type"].compare("edge") == 0)
                            {
                                if(j==0 || i==0 || i==x_size-1 || j==y_size-1) a_type = data.in.p_vec_int["Atom_type"][1];
                            }
                            if(data.in.param["Type"].compare("8point") == 0)
                            {
                                for(int index=0; index < 8; ++index)
                                {
                                    if( ((i-x[index])*(i-x[index]) + (j-y[index])*(j-y[index])) < rr*rr+0.1 )
                                    {
                                        a_type = data.in.p_vec_int["Atom_type"][1];
                                    }
                                }
                            }
                        }
                        if(data.in.p_vec_int["Atom_type"].size() == 5) // edge type, each edge different atom type
                        {
                            if(i == 0)                                   a_type = data.in.p_vec_int["Atom_type"][1];
                            if(i == x_size-1)                            a_type = data.in.p_vec_int["Atom_type"][2];
                            if(j == 0     && i != 0 && i != x_size-1)    a_type = data.in.p_vec_int["Atom_type"][3];
                            if(j == y_size-1 && i != 0 && i != x_size-1) a_type = data.in.p_vec_int["Atom_type"][4];
                        }
                    }



                    slab.push_back(Atom(N, pos, a_type, m_tag));
                    ++N;
                }
            }
        }


        if(data.in.p_float.contains("Radius"))
        {
            double radius = data.in.p_float["Radius"];
            double x,y;

            for(Atom& a : slab)
            {
                x = a.pos.x + (1.0/data.in.p_float["Scale"]) * data.in.p_tensor["Position_shift"].x;
                y = a.pos.y + (1.0/data.in.p_float["Scale"]) * data.in.p_tensor["Position_shift"].y;
                if(data.in.p_int.contains("Exclude_type") && data.in.p_int.contains("Exclude_mol_tag") && x*x + y*y < radius*radius)
                {
                    a.type = data.in.p_int["Exclude_type"];
                    a.mol_tag = data.in.p_int["Exclude_mol_tag"];
                }
            }
        }

        beads.insert(beads.end(), slab.begin(), slab.end());
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
        ss << "Slab_size: 10 10 2\n";
        ss << "Mol_tag: 1\n";
        ss << "Atom_type: 1\n";
        ss << "Atom_type: (3 options):\n";
        ss << "define single type for body: 1\n";
        ss << "define 2 types: 1->body 2->edges\n";
        ss << "define 5 types: 1->body 2, 3, 4, 5->for each edge\n";

        return ss.str();
    }

    void validate_inputs( Data& data )
    {
        data.in.p_tensor.validate_keyword("Position_shift", "0 0 0");
        data.in.p_float.validate_keyword("Scale", "1.0");
        data.in.p_int.validate_keyword("Mol_tag", "1");
        data.in.p_vec_int.validate_keyword("Atom_type", "1 2");
        data.in.p_vec_int.validate_keyword("Bond_type", "1 2");
        data.in.p_vec_int.validate_keyword("Slab_size", "10 10 2");

        if(data.in.p_vec_int["Slab_size"].size() != 3)
        {
            cerr << "Slab_size need to define 3 numbers. Slab_size: 10 10 2" << endl;
            exit(-1);
        }
    }

    void generate( Data& data )
    {
        validate_inputs(data);

        Atoms slab;
        int x_size = data.in.p_vec_int["Slab_size"][0];
        int y_size = data.in.p_vec_int["Slab_size"][1];
        int z_size = data.in.p_vec_int["Slab_size"][2];

        int N = 1; // offset handled by make persistent
        Tensor_xyz pos;
        int m_tag = data.in.p_int["Mol_tag"];
        int a_type = -1;

        for(int k=0; k<z_size; ++k)
        {
            for(int i=0; i<x_size; ++i)
            {
                for(int j=0; j<y_size; ++j)
                {
                    pos = Tensor_xyz(i, j, k);
                    a_type = data.in.p_vec_int["Atom_type"][0];

                    if(data.in.p_vec_int["Atom_type"].size() == 2 && k==0) // edges
                    {
                        if(j==0 || i==0 || i==x_size-1 || j==y_size-1) a_type = data.in.p_vec_int["Atom_type"][1];
                    }
                    if(data.in.p_vec_int["Atom_type"].size() == 5 && k==0)
                    {
                        if(i == 0)                                   a_type = data.in.p_vec_int["Atom_type"][1];
                        if(i == x_size-1)                            a_type = data.in.p_vec_int["Atom_type"][2];
                        if(j == 0     && i != 0 && i != x_size-1)    a_type = data.in.p_vec_int["Atom_type"][3];
                        if(j == y_size-1 && i != 0 && i != x_size-1) a_type = data.in.p_vec_int["Atom_type"][4];
                    }

                    slab.push_back(Atom(N, pos, a_type, m_tag));
                    ++N;
                }
            }
        }

        double radius = 0.0;
        double x,y;
        if(data.in.p_float.contains("Radius"))
        {
            radius = data.in.p_float["Radius"];
        }
        for(Atom& a : slab)
        {
            x = a.pos.x + (1.0/data.in.p_float["Scale"]) * data.in.p_tensor["Position_shift"].x;
            y = a.pos.y + (1.0/data.in.p_float["Scale"]) * data.in.p_tensor["Position_shift"].y;
            if(data.in.p_int.contains("Exclude_type") && data.in.p_int.contains("Exclude_mol_tag") && x*x + y*y < radius*radius)
            {
                a.type = data.in.p_int["Exclude_type"];
                a.mol_tag = data.in.p_int["Exclude_mol_tag"];
            }
        }

        gen_bonds(data, slab);

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

    void gen_bonds(Data& data, Atoms& slab)
    {
        int part_ID;
        int part_ID_i, part_ID_j, part_ID_k;
        int part_ID_ij, part_ID_ik, part_ID_jk;

        int s_x = data.in.p_vec_int["Slab_size"][0];
        int s_y = data.in.p_vec_int["Slab_size"][1];
        int s_z = data.in.p_vec_int["Slab_size"][2];
        vector<int> type_gen = data.in.p_vec_int["Atom_type"];
        vector<int> bond_type = data.in.p_vec_int["Bond_type"];

        for(int k=0; k<s_z; ++k)
        {
            for(int i=0; i<s_x; ++i)
            {
                for(int j=0; j<s_y; ++j)
                {
                    part_ID    = (k  )*s_x*s_y + (i  )*s_y + (j  );
                    part_ID_k  = (k+1)*s_x*s_y + (i  )*s_y + (j  );
                    part_ID_i  = (k  )*s_x*s_y + (i+1)*s_y + (j  );
                    part_ID_j  = (k  )*s_x*s_y + (i  )*s_y + (j+1);
                    part_ID_ij = (k  )*s_x*s_y + (i+1)*s_y + (j+1);
                    part_ID_ik = (k+1)*s_x*s_y + (i+1)*s_y + (j  );
                    part_ID_jk = (k+1)*s_x*s_y + (i  )*s_y + (j+1);

                    if(i < s_x-1 && is_bond_type(slab[part_ID].type, slab[part_ID_i].type, type_gen) ) // bond in i dir
                    {
                        bonds.push_back(    Bond(bonds.size()+1, bond_type[0], slab[part_ID_i].N, slab[part_ID].N)    );
                    }
                    if(j < s_y-1 && is_bond_type(slab[part_ID].type, slab[part_ID_j].type, type_gen) ) // bond in j dir
                    {
                        bonds.push_back(    Bond(bonds.size()+1, bond_type[0], slab[part_ID_j].N, slab[part_ID].N)    );
                    }
                    if(k < s_z-1 && is_bond_type(slab[part_ID].type, slab[part_ID_k].type, type_gen) ) // bond in k dir
                    {
                        bonds.push_back(    Bond(bonds.size()+1, bond_type[0], slab[part_ID_k].N, slab[part_ID].N)    );
                    }

                    if(bond_type.size() >= 2)
                    {
                        if(i < s_x-1 && j < s_y-1 && is_bond_type(slab[part_ID].type, slab[part_ID_ij].type, type_gen) ) // bond in ij plane
                        {
                            bonds.push_back(    Bond(bonds.size()+1, bond_type[1], slab[part_ID_ij].N, slab[part_ID].N)    );
                        }
                        if(i < s_x-1 && k < s_z-1 && is_bond_type(slab[part_ID].type, slab[part_ID_ik].type, type_gen) ) // bond in ik plane
                        {
                            bonds.push_back(    Bond(bonds.size()+1, bond_type[1], slab[part_ID_ik].N, slab[part_ID].N)    );
                        }
                        if(j < s_y-1 && k < s_z-1 && is_bond_type(slab[part_ID].type, slab[part_ID_jk].type, type_gen) ) // bond in jk plane
                        {
                            bonds.push_back(    Bond(bonds.size()+1, bond_type[1], slab[part_ID_jk].N, slab[part_ID].N)    );
                        }
                    }
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
