#ifndef VIRUS_PSEUDOT3_H
#define VIRUS_PSEUDOT3_H

#include "system_base.h"
#include "atom.h"

class T3_Penta
{
public:
    int from[5];
    int to[5];

    Atoms coms;

    friend std::ostream& operator<<(std::ostream& os, T3_Penta& pp)
    {
        for(int i=0; i<pp.coms.size(); ++i)
        {
            os << "Com " << i << " = " << pp.coms[i] << " index from " << pp.from[i] << " to " << pp.to[i] << endl;
        }
        return os;
    }
};

class Virus_pseudoT3 : public System_Base
{
public:
    inline static const string keyword = "Virus_T=3p";
    const string name = "Virus_T=3p";

    Virus_pseudoT3() : System_Base("Virus_T=3p") {}

    void generate( Data& data )
    {
        if(data.in.id == 2)
        {
            set_protomer(data.beads[0]); // TODO: make it work with input file
            set_capsid(data.beads[data.id_map[data.in.id]]); // TODO: make it work with input file
            calc_protomer_coms();
            calc_pentamers();
            calc_symmetry_axes();

            align_2fold_on_x();
        }
    }

    string help()
    {
        stringstream ss;

        ss << "Load 2 files, first containing protomer, then full capsid" << endl;
        ss << "System_type: Protomer\n";
        ss << "Output_type: pdb # other keywords: xyz pdb lammps_full\n";

        ss << "System_type: Virus_T=3*\n";
        ss << "Output_type: pdb # other keywords: xyz pdb lammps_full\n";

        return ss.str();
    }

    Atoms protomer;
    Atoms capsid;
    Atoms proto_coms;

    vector<T3_Penta> penta;

    Atoms sym_axis_2;
    Atoms sym_axis_3;
    Atoms sym_axis_5;

    void printXYZ()
    {
        cout << sym_axis_2.size() + sym_axis_3.size() + sym_axis_5.size() << "\nparticle\n";
        for (Atom& a : sym_axis_2)
        {
            cout << "C" << " " << a.pos << "\n";
        }
        for (Atom& a : sym_axis_3)
        {
            cout << "O" << " " << a.pos << "\n";
        }
        for (Atom& a : sym_axis_5)
        {
            cout << "N" << " " << a.pos << "\n";
        }
    }

    void set_protomer(Atoms& a)
    {
        protomer = a;
    }

    void set_capsid(Atoms& a)
    {
        capsid = a;
    }

    void align_2fold_on_x()
    {
        //cout << "exit in Data::align_2fold_on_x" << endl;
        //exit(1);
    }

    void calc_protomer_coms()
    {
        if(protomer.empty() || capsid.empty())
        {
            cerr << "Virus_pseudoT3::calc_protomer_coms executed without protomer or capsid set" << endl;
            exit(1);
        }
        proto_coms = get_protomer_coms();
    }

    void calc_pentamers()
    {
        bool exists = false;
        int id;

        //
        // store identified pentamers
        //
        for(Atom& com : proto_coms)
        {
            exists = false;
            for(T3_Penta& p : penta)
            {
                for(Atom& cm : p.coms)
                {
                    if(cm.isAproxSame(com))
                    {
                        exists = true;
                    }
                }
            }
            if(!exists)
            {
                penta.push_back(get_pentamer(com));
            }
        }

        //
        // Calc protomer index ranges of pentamers
        //
        for(T3_Penta& p : penta)
        {

            for(int i=0; i<p.coms.size(); ++i)
            {
                id = get_protomer_id( p.coms[i] );
                p.from[i] = id*protomer.size();
                p.to[i] = (id+1)*protomer.size();
            }
        }

        //
        // Test
        //
        if(penta.size() != 12)
        {
            cerr << "Virus_pseudoT3::calc_pentamers identified " << penta.size() << " pentamers" << endl;
            exit(1);
        }
    }


    void calc_symmetry_axes()
    {
        get_5_axis();
        get_2_axis();
        get_3_axis();
    }

private:

    Atoms get_protomer_coms()
    {
        Atoms com;
        Atom cm;

        for(int i=0; i<60; ++i)
        {
            cm = get_protomer_com(i);
            com.push_back(cm);
        }

        return com;
    }

    int get_protomer_id(Atom com)
    {
        for(int i=0; i<60; ++i)
        {
            if( com.isAproxSame(get_protomer_com(i)) )
            {
                return i;
            }
        }

        cerr << "Virus_pseudoT3::get_protomer_id failed to get protomer id" << endl;
        exit(-1);
    }

    Atom get_protomer_com(int i)
    {
        return capsid.center_of_mass( i*protomer.size(), (i+1)*protomer.size() );
    }

    void get_5_axis()
    {
        Atoms sa5;
        for(T3_Penta& p : penta)
        {
            sa5.push_back( p.coms.center_of_mass() );
        }

        sa5.scale( proto_coms[0].size() / sa5[0].size() );

        if(sa5.size() != 12)
        {
            cerr << "Virus_pseudoT3::get_5_axis identified " << sa5.size() << " symmetry axes" << endl;
            exit(1);
        }

        sym_axis_5 = sa5;
    }

    void get_2_axis()
    {
        if(sym_axis_5.empty())
        {
            exit(1);
        }

        bool exists = false;
        double min_d = sym_axis_5.min_dist();
        Atoms sa2;
        Atoms cm;

        //
        // 2-fold symmetry axis == COM of closest pair of 5-fold symmetry axis
        //
        for(Atom& i : sym_axis_5)
        {
            for(Atom& j : sym_axis_5)
            {
                if(i != j && i.dist(j) < min_d*1.05)
                {
                    cm.push_back(i);
                    cm.push_back(j);

                    exists = false;
                    for(auto& a : sa2)
                    {
                        if(a.isAproxSame(cm.center_of_mass()))
                            exists = true;
                    }

                    if(!exists)
                        sa2.push_back( cm.center_of_mass() );

                    cm.clear();
                }
            }
        }

        sa2.scale( proto_coms[0].size() / sa2[0].size() );

        if(sa2.size() != 30)
        {
            cerr << "Virus_pseudoT3::get_2_axis identified " << sa2.size() << " symmetry axes" << endl;
            exit(1);
        }

        sym_axis_2 = sa2;
    }

    void get_3_axis()
    {
        bool exists = false;
        double min_d = sym_axis_5.min_dist();
        Atoms sa3;
        Atoms cm;

        for(auto& i : sym_axis_5)
        {
            for(auto& j : sym_axis_5)
            {
                for(auto& k : sym_axis_5)
                {
                    if( i != j && i.dist(j) < min_d*1.05 &&
                        i != k && i.dist(k) < min_d*1.05 &&
                        j != k && j.dist(k) < min_d*1.05 )
                    {
                        cm.push_back(i);
                        cm.push_back(j);
                        cm.push_back(k);

                        exists = false;
                        for(auto& a : sa3)
                        {
                            if(a.isAproxSame(cm.center_of_mass()))
                                exists = true;
                        }

                        if(!exists)
                        {
                            sa3.push_back( cm.center_of_mass() );
                        }

                        cm.clear();
                    }
                }
            }
        }

        sa3.scale( proto_coms[0].size() / sa3[0].size() );

        if(sa3.size() != 20)
        {
            cerr << "Virus_pseudoT3::get_3_axis identified " << sa3.size() << endl;
            exit(1);
        }

        sym_axis_3 = sa3;
    }

    T3_Penta get_pentamer(Atom origin)
    {
        T3_Penta penta;
        Atoms sel;

        // Protomer neighbors in a pentamer should have shortest distance
        double min_d = proto_coms.min_dist();

        sel = proto_coms.within(origin, 1.01*min_d);
        for(Atom& item : sel)
        {
            sel.join( proto_coms.within(item, 1.01*min_d) );
        }

        if(sel.size() != 5)
        {
            cerr << "Virus_pseudoT3::get_pentamer() failed to get a pentamer" << endl;
            exit(1);
        }

        penta.coms = sel;
        return penta;
    }
};

#endif // VIRUS_PSEUDOT3_H
