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

using namespace std;

/**
 * @brief ran - 0.0 to 1.0
 * @return
 */
double ran() {
    return unif(rng);
}














class Quat{
public:
    myFloat w,x,y,z;

    Quat(myFloat w,myFloat x, myFloat y, myFloat z): w(w), x(x), y(y), z(z) {}
};




class Tensor_xyz
{
public:
	Tensor_xyz(){}
	Tensor_xyz(myFloat x, myFloat y, myFloat z): x(x), y(y), z(z) {}

	myFloat x=0.0, y=0.0, z=0.0;


    double size() const {
        return sqrt(x*x + y*y + z*z);
    }

    /**
     * @brief Normalise a vector to have unit length.
     * For speed during heavy use, it is not checked that the supplied vector has non-zero length.
     */
    inline void normalise()
    {
        myFloat tot = size();
        if (tot !=0.0) {
            tot = 1.0 / tot;
            x *= tot;
            y *= tot;
            z *= tot;
        }
    }

    inline void randomUnitSphere()
    {
        myFloat a, b, xi1, xi2;

        do {
            xi1 = 1.0 - 2.0 * ( ran() );
            xi2 = 1.0 - 2.0 * ( ran() );

            a = xi1*xi1 + xi2*xi2;
        } while (a > 1.0);

        b = 2.0 * sqrt(1.0 - a);

        x = xi1 * b;
        y = xi2 * b;
        z = 1.0 - 2.0*a;
    }

    bool isAproxSame(const Tensor_xyz& o, myFloat approx = 0.000001) const {

        if(!(x < o.x+approx && x > o.x - approx))
            return false;
        if(!(y < o.y+approx && y > o.y - approx))
            return false;
        if(!(z < o.z+approx && z > o.z - approx))
            return false;

        return true;
    }

    double dot(const Tensor_xyz& o) const {
        return this->x*o.x + this->y*o.y + this->z*o.z;
    }

    inline Tensor_xyz cross(const Tensor_xyz& B) const {
        return Tensor_xyz(y*B.z - z*B.y, -x*B.z + z*B.x, x*B.y - y*B.x);
    }

    double dist(const Tensor_xyz& o) const {
        return sqrt( (this->x-o.x)*(this->x-o.x) + (this->y-o.y)*(this->y-o.y) + (this->z-o.z)*(this->z-o.z) );
    }

    inline double distSQ(const Tensor_xyz& o) const {
        return (this->x-o.x)*(this->x-o.x) + (this->y-o.y)*(this->y-o.y) + (this->z-o.z)*(this->z-o.z);
    }

    inline bool checkBox(myFloat size) {
        return x<size && x>-size && y<size && y>-size && z<size && z>-size;
    }

    inline void rotate(Tensor_xyz& axis, myFloat angle) {
        angle*=0.5;
        double cosAngle = cos(angle);
        double sinAngle = sin(angle);
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz;
        double qw = cosAngle, qx = (axis.x * sinAngle), qy = (axis.y * sinAngle), qz = (axis.z * sinAngle);

        /*    t1 = quat.w * quat.w; */
        t2 =  qw * qx;
        t3 =  qw * qy;
        t4 =  qw * qz;
        t5 = -qx * qx;
        t6 =  qx * qy;
        t7 =  qx * qz;
        t8 = -qy * qy;
        t9 =  qy * qz;
        t10 = -qz * qz;

        newx = 2.0 * ( (t8+t10) * x + (t6-t4)  * y + (t3+t7) * z ) + x;
        newy = 2.0 * ( (t4+t6)  * x + (t5+t10) * y + (t9-t2) * z ) + y;
        newz = 2.0 * ( (t7-t3)  * x + (t2+t9)  * y + (t5+t8) * z ) + z;

        x = newx;
        y = newy;
        z = newz;
    }

    inline void rotate(Quat& q) {
        double t2,t3,t4,t5,t6,t7,t8,t9,t10,newx,newy,newz;

        /*    t1 = quat.w * quat.w; */
        t2 =  q.w * q.x;
        t3 =  q.w * q.y;
        t4 =  q.w * q.z;
        t5 = -q.x * q.x;
        t6 =  q.x * q.y;
        t7 =  q.x * q.z;
        t8 = -q.y * q.y;
        t9 =  q.y * q.z;
        t10 = -q.z * q.z;

        newx = 2.0 * ( (t8+t10) * x + (t6-t4)  * y + (t3+t7) * z ) + x;
        newy = 2.0 * ( (t4+t6)  * x + (t5+t10) * y + (t9-t2) * z ) + y;
        newz = 2.0 * ( (t7-t3)  * x + (t2+t9)  * y + (t5+t8) * z ) + z;

        x = newx;
        y = newy;
        z = newz;
    }

    //
    //
    // Operators
    //
    //
    bool operator==(const Tensor_xyz& o) const
    {
        if(x != o.x) return false;
        if(y != o.y) return false;
        if(z != o.z) return false;
        /*const double aprox = 0.0000000000001;
        if(x < o.x+aprox && x > o.x - aprox)
            return false;
        if(y < o.y+aprox && y > o.y - aprox)
            return false;
        if(z < o.z+aprox && z > o.z - aprox)
            return false;*/

        return true;
    }

    void operator*=(myFloat a) {
        x*=a;
        y*=a;
        z*=a;
    }

    void operator/=(myFloat a) {
        x/=a;
        y/=a;
        z/=a;
    }

    void operator+=(const Tensor_xyz& o) {
        x += o.x;
        y += o.y;
        z += o.z;
    }

    Tensor_xyz operator*(const myFloat a) const {
        return Tensor_xyz(x*a,y*a,z*a);
    }

    Tensor_xyz operator/(const myFloat a) const {
        return Tensor_xyz(x/a,y/a,z/a);
    }

    Tensor_xyz operator+(const Tensor_xyz& o) const {
        return Tensor_xyz(x+o.x, y+o.y, z+o.z);
    }

    Tensor_xyz operator-(const Tensor_xyz& o) const {
        return Tensor_xyz(x-o.x, y-o.y, z-o.z);
    }

    friend std::ostream& operator<<(std::ostream& os, const Tensor_xyz& vec) {
      os << vec.x << " " << vec.y << " " << vec.z;
      return os;
    }
};




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
        return ( lj[type1].sigma + lj[type2].sigma );
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
