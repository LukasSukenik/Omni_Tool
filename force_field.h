#ifndef FORCE_FIELD_H
#define FORCE_FIELD_H

#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <cmath>

#include "rng.h"
#include "types.h"
#include "angle.h"
#include "bond.h"
#include "tensor_xyz.h"

using namespace std;


























/**
 * @brief The BeadParam class
 * Bead LJ parameters
 */
class LJ {
public:
	LJ(){}
	LJ(int type, double epsilon, double sigma, double cutoff) : type(type), epsilon(epsilon), sigma(sigma), cutoff(cutoff) {}
	LJ(int type1, int type2, double epsilon, double sigma, double cutoff) : type1(type1), type2(type2), epsilon(epsilon), sigma(sigma), cutoff(cutoff) {}

    int type=-1;
    int type1=-1;
    int type2=-1;
    double epsilon=1.0;
    double sigma=1.0;
    double cutoff=25.0;

    double energy(double r)
    {
    	if(r > cutoff)
    		return 0.0;
    	return 4*epsilon * ( pow(sigma/r, 12) - pow(sigma/r, 6) );
    }

    double force(double r)
    {
    	if(r > cutoff)
    		return 0.0;
    	return 24*epsilon * ( 2*pow(sigma/r, 13) - pow(sigma/r, 7) );
    }

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
    int type1=-1;
    int type2=-1;
    double epsilon=1.0;
    double start_dis=1.0;
    double range=1.0;
    double pi = 3.141592653589793238462643383279502884197;

    //
    // Energy: epsilon * cos^2 ( pi/2 * ( r - start_dis) / range )
    //
    // Force: - epsilon*pi/range * sin(arg)*cos(arg)
    // arg = 0.5*pi * (r - start_dis) / range
    //

    double energy(double r)
    {
    	if(r < start_dis || r > start_dis + range)
    	{
    	  	return 0.0;
    	}

    	return - epsilon * pow( cos(0.5*pi * (r - start_dis) / range), 2);
    }

    double force(double r)
    {
    	if(r < start_dis || r > start_dis + range)
    	{
    	  	return 0.0;
    	}

    	double argument = 0.5*pi*(r-start_dis)/range;
    	double sinn, coss;
    	sincos(argument, &sinn, &coss);

    	return ( epsilon*pi/range ) * coss*sinn;
    }

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
    Force_Field() {}

    vector<int> types;
    map<int, LJ> lj;
    map<int,CosSQ> cos;

    bool empty()
    {
        return types.empty() && lj.empty() && cos.empty();
    }

    double energy(double r, int type1, int type2)
    {
    	//return cos[type1].energy(r);
    	return lj[type1].energy(r) + cos[type1].energy(r);
    }

    double force(double r, int type1, int type2)
    {
    	//return cos[type1].force(r);
    	return lj[type1].force(r) + cos[type1].force(r);
    }

    double get_cutoff(int type1, int type2)
    {
        //cerr << type1 << ":" << type2 << " = " << 0.5 * ( lj[type1].sigma + lj[type2].sigma ) << endl;
        return 0.5 * ( lj[type1].sigma + lj[type2].sigma );
    }

    friend std::ostream& operator<<(std::ostream& os, Force_Field& ff)
    {
        os << "\n";
        for(auto& [key, val] : ff.lj)
        {
            os << "LJ" << key << ": " << val << "\n";
        }
        for(auto& [key, val] : ff.cos)
        {
            os << "CosSQ" << key << ": " << val << "\n";
        }
        return os;
    }

};

#endif // FORCE_FIELD_H
