#ifndef FORCE_FIELD_H
#define FORCE_FIELD_H

#include <vector>
#include <string>

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
    int calcBondTypes() const
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
};

#endif // FORCE_FIELD_H
