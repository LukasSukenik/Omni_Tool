#ifndef VECTOR_H
#define VECTOR_H

#include <string>

#include "rng.h"
#include "force_field.h"


double ran() {
    return unif(rng);
}




class Quat{
public:
    myFloat w,x,y,z;

    Quat(myFloat w,myFloat x, myFloat y, myFloat z): w(w), x(x), y(y), z(z) {}
};




class Atom{
public:
    //
    // Lammps define Atom/Particle
    //
    myFloat x=0,y=0,z=0;
    myFloat vx=-1.001,vy=-1.001,vz=-1.001; // velocities
    double q=0;
    int nx=0,ny=0,nz=0;

    /*
     * From pdb: Atom name < Residue name < Chain identifier
     */

    int N=1; // PDB: Atom serial number
    int type=0; // PBD: Atom name
    std::string resname = "   "; // PDB: Residue name
    int mol_tag = 0; // PDB: Chain identifier and Residue sequence number


    //
    // not used for output
    //
    int to_type=-1; // used in dodecahedron to change the type after bonds are generated - design mistake
    double temp_dist = -1.0;

    //
    // Constructors
    //
    Atom() : to_type(-1) {}
    Atom(myFloat x, myFloat y, myFloat z, int type=0): x(x), y(y), z(z), type(type), to_type(-1) {}
    Atom(myFloat x, myFloat y, myFloat z, myFloat vx, myFloat vy, myFloat vz, int type=0): x(x), y(y), z(z), vx(vx), vy(vy), vz(vz), type(type), to_type(-1) {}
    Atom(myFloat x, myFloat y, myFloat z, int type, int mol_tag): x(x), y(y), z(z), type(type), mol_tag(mol_tag), to_type(-1) {}

    bool operator==(const Atom& o) const {
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

    bool operator!=(const Atom& o) const {
        return !(*this == o);
    }

    /**
     * @brief normalise Normalise a vector to have unit length.  For speed during heavy use, it is
       not checked that the supplied vector has non-zero length.
     */
    inline void normalise() {
        myFloat tot = size();
        if (tot !=0.0) {
            tot = 1.0 / tot;
            x *= tot;
            y *= tot;
            z *= tot;
        }
    }

    inline void randomUnitSphere() {
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

    bool isAproxSame(const Atom& o, myFloat approx = 0.000001) const {

        if(!(x < o.x+approx && x > o.x - approx))
            return false;
        if(!(y < o.y+approx && y > o.y - approx))
            return false;
        if(!(z < o.z+approx && z > o.z - approx))
            return false;

        return true;
    }

    double size() const {
        return sqrt(x*x + y*y + z*z);
    }

    Atom operator*(const myFloat a) const {
        return Atom(x*a,y*a,z*a, type, mol_tag);
    }

    Atom operator/(const myFloat a) const {
        return Atom(x/a,y/a,z/a, type, mol_tag);
    }

    void operator*=(myFloat a) {
        x*=a;
        y*=a;
        z*=a;
    }

    double dot(const Atom& o) const {
        return this->x*o.x + this->y*o.y + this->z*o.z;
    }

    Atom operator+(const Atom& o) const {
        return Atom(x+o.x, y+o.y, z+o.z, type, mol_tag);
    }

    void operator+=(const Atom& o) {
        x += o.x;
        y += o.y;
        z += o.z;
    }

    Atom operator-(const Atom& o) const {
        return Atom(x-o.x, y-o.y, z-o.z, type);
    }

    double dist(const Atom& o) const {
        return sqrt( (this->x-o.x)*(this->x-o.x) + (this->y-o.y)*(this->y-o.y) + (this->z-o.z)*(this->z-o.z) );
    }

    inline double distSQ(const Atom& o) const {
        return (this->x-o.x)*(this->x-o.x) + (this->y-o.y)*(this->y-o.y) + (this->z-o.z)*(this->z-o.z);
    }

    inline bool checkBox(myFloat size) {
        return x<size && x>-size && y<size && y>-size && z<size && z>-size;
    }

    inline Atom cross(const Atom& B) const {
        return Atom(this->y*B.z - this->z*B.y, -this->x*B.z + this->z*B.x, this->x*B.y - this->y*B.x);
    }

    bool isNeighbor(const Atom& o, myFloat len = 1.0, myFloat margin = 0.00000001) const
    {
        Atom vec(o.x - x, o.y - y, o.z - z);
        double dist = vec.x*vec.x + vec.y*vec.y + vec.z*vec.z;
        return ( dist < len+margin && dist > len-margin );
    }

    inline void rotate(Atom& axis, myFloat angle) {
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


    friend std::ostream& operator<<(std::ostream& os, const Atom& vec) {
      os << vec.x << " " << vec.y << " " << vec.z;
      return os;
    }
};




bool isAproxSame(const myFloat& a, const myFloat& b, myFloat approx = 0.000001) {
    return (a < b+approx && a > b - approx);
}




bool myfunction (Atom i,Atom j) { return (i.size()<j.size()); }




class Atoms : public vector< Atom >
{
public:
	//
	//
	// Methods that change something in the container
	//
	//

	void set_mol_tag(int mtag)
	{
	    for(Atom& item : (*this))
	        item.mol_tag = mtag;
	}

    void move(Atom move)
    {
        for(Atom& item : (*this))
            item += move;
    }

    void scale(double scale)
    {
        for(Atom& item : (*this))
            item *= scale;
    }

    void project_to_unit_sphere()
    {
    	Atom zero = Atom(0.0,0.0,0.0);

    	for(Atom& a : (*this))
    	{
    		a *= 1.0/a.dist(zero);
    	}
    }

    //
    //
    // Methods const
    //
    //

    int count_Atom_Type( int atype) const
    {
        int count = 0;
        for(const Atom& item : (*this))
        {
            if( item.type == atype)
            {
                ++count;
            }
        }
        return count;
    }

    int count_Mol_tag(int mTag) const
    {
        int count = 0;
        for(const Atom& item : (*this))
        {
            if( item.mol_tag == mTag)
            {
                ++count;
            }
        }
        return count;
    }

    vector<int> get_Atom_Types() const
    {
        vector<int> atom_types;
        bool exist = false;

        atom_types.push_back((*this)[0].type);
        for(const Atom& a : (*this))
        {
            exist = false;
            for(int j=0; j<atom_types.size(); ++j)
            {
                if(atom_types[j] == a.type)
                    exist = true;
            }
            if(!exist)
                atom_types.push_back( a.type );
        }
        return atom_types;
    }

    int getMaxMolTag()
    {
        int max=0;
        for(auto& a : (*this))
        {
            if(a.mol_tag > max)
            {
                max = a.mol_tag;
            }
        }
        return max;
    }

    bool is_overlap(Atoms other, Force_Field ff) const
    {
    	for(const Atom& a : (*this))
    	{
    		for(const Atom& o : other)
    		{
                if(a.distSQ(o) < ff.get_cutoff(a.type, o.type)*ff.get_cutoff(a.type, o.type))
    				return true;
    		}
    	}
    	return false;
    }

    //
    //
    // static methods
    //
    //

    static Atom clusterCM(vector<Atom>& cluster) {

        Atom cluscm(0.0, 0.0, 0.0);

        for(unsigned int i=0; i<cluster.size(); i++) {
            cluscm.x += cluster[i].x;
            cluscm.y += cluster[i].y;
            cluscm.z += cluster[i].z;
        }

        cluscm.x /= cluster.size();
        cluscm.y /= cluster.size();
        cluscm.z /= cluster.size();

        return cluscm;
    }

    static void clusterRotate(vector<Atom>& cluster, double angle, Atom axis) {
        double vc,vs;

        Atom cluscm = clusterCM(cluster);

        axis.normalise();

        vc = cos( angle );
        vs = sqrt(1.0 - vc*vc);

        Quat newquat(vc, axis.x*vs, axis.y*vs, axis.z*vs);

        //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

        //shift position to geometrical center
        for(unsigned int i=0; i<cluster.size(); ++i) {
            //shift position to geometrical center
            cluster[i].x -= cluscm.x;
            cluster[i].y -= cluscm.y;
            cluster[i].z -= cluscm.z;
            //do rotation
            cluster[i].rotate(newquat);
            //shift positions back
            cluster[i].x += cluscm.x;
            cluster[i].y += cluscm.y;
            cluster[i].z += cluscm.z;
        }
    }

    static void clusterRotate_random(vector<Atom>& cluster, double max_angle) {
        double vc,vs;
        Atom newaxis;

        Atom cluscm = clusterCM(cluster);

        // create rotation quaternion
        newaxis.randomUnitSphere(); // random axes for rotation
        vc = cos(max_angle * ( 0.0001 * (rand()%10000) ) );
        if (( 0.0001 * (rand()%10000) ) <0.5) vs = sqrt(1.0 - vc*vc);
        else vs = -sqrt(1.0 - vc*vc); // randomly choose orientation of direction of rotation clockwise or counterclockwise

        Quat newquat(vc, newaxis.x*vs, newaxis.y*vs, newaxis.z*vs);

        //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

        //shift position to geometrical center
        for(unsigned int i=0; i<cluster.size(); ++i) {
            //shift position to geometrical center
            cluster[i].x -= cluscm.x;
            cluster[i].y -= cluscm.y;
            cluster[i].z -= cluscm.z;
            //do rotation
            cluster[i].rotate(newquat);
            //shift positions back
            cluster[i].x += cluscm.x;
            cluster[i].y += cluscm.y;
            cluster[i].z += cluscm.z;
        }
    }
};

#endif // VECTOR_H
