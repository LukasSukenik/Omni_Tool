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
    }

    void generate( Data& data )
    {
        validate_inputs(data);

        int row = sqrt(data.in.p_int["Number_of_beads"]);
        int N = 1;

        for(int i=0; i<row; ++i) {
            for(int j=0; j<row; ++j) {
                for(int k=0; k<1; ++k) {
                    if( j==0 || i == 0 || i == row-1 || j == row-1 )
                        beads.push_back(Atom(N, Tensor_xyz(i, k, j), data.in.p_vec_int["Atom_type"][1], data.in.p_int["Mol_tag"]));
                    else
                        beads.push_back(Atom(N, Tensor_xyz(i, k, j), data.in.p_vec_int["Atom_type"][0], data.in.p_int["Mol_tag"]));
                    ++N;
                }
            }
        }

        gen_bonds(row);
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

    void gen_bonds(int row) {
        int actual;
        for(int i=0; i<row; ++i) {
            for(int j=0; j<row; ++j) {
                actual = i*row + j;
                if(j > 0)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i)*row + j-1, actual, beads[(i)*row + j-1].dist(beads[actual]) ) );
                if(i > 0)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i-1)*row + j, actual, beads[(i-1)*row + j].dist(beads[actual]) ) );
                if(j < row-1)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i)*row + j+1, actual, beads[(i)*row + j+1].dist(beads[actual]) ) );
                if(i < row-1)
                    bonds.push_back( Bond(bonds.size(), bonds.size(), (i+1)*row + j, actual, beads[(i+1)*row + j].dist(beads[actual]) ) );
            }
        }

        // Correct for offset -> lammps starts at 1 not 0
        for(int i=0; i<bonds.size(); ++i) {
            ++bonds[i].at1;
            ++bonds[i].at2;
        }
    }
};

#endif // SLAB_H
