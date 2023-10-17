#ifndef FORCE_FIELD_H
#define FORCE_FIELD_H

#include <sstream>
#include <vector>
#include <string>
#include <iostream>

using namespace std;

class Bond {
public:
    Bond() : N(-1), type(-1), at1(-1), at2(-1) {}
    Bond(int N, int type, int at1, int at2) : N(N), type(type), at1(at1), at2(at2) {}
    Bond(int N, int type, int at1, int at2, bool typelock) : N(N), type(type), at1(at1), at2(at2), typelock(typelock) {}
    Bond(int N, int type, int at1, int at2, double r0) : N(N), type(type), at1(at1), at2(at2), r0(r0), coeff_name("coeff_bond") {}
    Bond(int N, int type, int at1, int at2, double r0, string coeff_name) : N(N), type(type), at1(at1), at2(at2), r0(r0), coeff_name(coeff_name) {}
    int N, type, at1, at2;
    bool typelock = false;

    double r0; /// equilibrium distance
    double K; /// energy/distance^2 -> 1/2 included

    string coeff_name;
};

bool sort_Bond_by_type(const Bond& i, const Bond& j) {
    return ( i.type < j.type ) || ( (i.at1 < j.at1) && (i.type == j.type) ) ;
}

class Bonds : public vector< Bond >
{
public:
    int calc_Bond_Types() const
    {
        vector<int> bond_types;
        bool exist = false;

        if(!empty())
        {
            bond_types.push_back((*this)[0].type);
            for(auto& bond : (*this)) {

                exist = false;
                for(auto btype : bond_types) {
                    if(btype == bond.type) {
                        exist = true;
                    }
                }
                if(!exist) {
                    bond_types.push_back( bond.type );
                }
            }
        }

        return bond_types.size();
    }

    void removeDuplicateBond()
    {
        int count = 0;
        for(auto& a : (*this))
        {
            for(auto& b : (*this))
            {
                if(a.at1 == b.at1 && a.at2 == b.at2 && a.N != b.N)
                {
                    cerr << "Duplicate bond " << a.N << "==" << b.N << endl;
                    ++count;
                }
            }
        }
        cerr << count << endl;
    }
};




class Angle {
public:
    Angle() : N(-1), type(-1), at1(-1), at2(-1), at3(-1) {}
    Angle(int N, int type, int at1, int at2, int at3) : N(N), type(type), at1(at1), at2(at2), at3(at3) {}
    int N, type, at1, at2, at3;
};




/**
 * @brief The BeadParam class
 * Bead LJ parameters
 */
class LJ {
public:
	LJ(){}
	LJ(int type, double epsilon, double sigma, double cutoff) : type(type), epsilon(epsilon), sigma(sigma), cutoff(cutoff) {}

    int type=-1;
    double epsilon=1.0;
    double sigma=1.0;
    double cutoff=25.0;

    friend std::ostream& operator<<(std::ostream& os, const LJ& lj)
    {
        os << lj.type << " " << lj.epsilon << " " << lj.sigma << " " << lj.cutoff;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, LJ& lj)
    {
        is >> lj.type >> lj.epsilon >> lj.sigma >> lj.cutoff;
        return is;
    }
};

class CosSQ {
public:
	CosSQ(){}
    CosSQ(int type, double epsilon, double start_dis, double range) : type(type), epsilon(epsilon), start_dis(start_dis), range(range) {}
	CosSQ(int type1, int type2, double epsilon, double start_dis, double range) : type1(type1), type2(type2), epsilon(epsilon), start_dis(start_dis), range(range){}

    int type=-1;
    int type1;
    int type2;
    double epsilon=1.0;
    double start_dis=1.0;
    double range=1.0;

    friend std::ostream& operator<<(std::ostream& os, const CosSQ& cos)
    {
        os << cos.type << " " << cos.epsilon << " " << cos.start_dis << " " << cos.range;
        return os;
    }

    friend std::istream& operator>>(std::istream& is, CosSQ& cos)
    {
        is >> cos.type >> cos.epsilon >> cos.start_dis >> cos.range;
        return is;
    }
};




class Force_Field
{
public:
    Force_Field() {
        lj.push_back(LJ(0,0,0,0)); // Types for lammps start at 1, but in c++ array starts at 0, so we fill the 0 position
        cos.push_back(CosSQ(0,0,0,0));
    }

    vector<int> types;
	vector<LJ> lj;
	vector<CosSQ> cos;

    double get_cutoff(int type1, int type2)
    {
        //cerr << type1 << ":" << type2 << " = " << 0.5 * ( lj[type1].sigma + lj[type2].sigma ) << endl;
        return ( lj[type1].sigma + lj[type2].sigma );
    }

    friend std::ostream& operator<<(std::ostream& os, const Force_Field& ff)
    {
        os << "\n";
        for(int i = 1; i < ff.lj.size(); ++i)
        {
            os << "LJ" << i << ": " << ff.lj[i] << "\n";
        }
        for(int i = 1; i < ff.cos.size(); ++i)
        {
            os << "CosSQ" << i << ": " << ff.cos[i] << "\n";
        }
        return os;
    }
};

#endif // FORCE_FIELD_H
