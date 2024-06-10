#ifndef VECTOR_H
#define VECTOR_H

#include <string>

#include "rng.h"
#include "force_field.h"







class Atom{
public:
    //
    // Lammps define Atom/Particle
    //
    //myFloat x=0,y=0,z=0;
    //myFloat vx=-1.001,vy=-1.001,vz=-1.001; // velocities

    Tensor_xyz pos = Tensor_xyz(0.0, 0.0, 0.0);
    Tensor_xyz vel = Tensor_xyz(1.001, 1.001, 1.001);

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
    Atom(Tensor_xyz pos, int type=0): pos(pos), type(type), to_type(-1) {}
    Atom(myFloat x, myFloat y, myFloat z, int type=0): pos(x,y,z), type(type), to_type(-1) {}

    Atom(myFloat x, myFloat y, myFloat z, myFloat vx, myFloat vy, myFloat vz, int type=0): pos(x,y,z), vel(vx,vy,vz), type(type), to_type(-1) {}

    Atom(myFloat x, myFloat y, myFloat z, int type, int mol_tag): pos(x,y,z), type(type), mol_tag(mol_tag), to_type(-1) {}
    Atom(Tensor_xyz pos, int type, int mol_tag): pos(pos), type(type), mol_tag(mol_tag), to_type(-1) {}

    bool operator==(const Atom& o) const
    {
    	return (this->pos == o.pos);
    }

    bool operator!=(const Atom& o) const {
        return !(*this == o);
    }

    void operator*=(myFloat a)
    {
    	pos *= a;
    }

    void operator/=(myFloat a)
    {
    	pos /= a;
    }

    void operator+=(const Atom& o)
    {
    	pos += o.pos;
    }

    Atom operator*(const myFloat a) const {
        return Atom(pos*a, type, mol_tag);
    }

    Atom operator/(const myFloat a) const {
        return Atom(pos/a, type, mol_tag);
    }

    Atom operator+(const Atom& o) const {
        return Atom(pos+o.pos, type, mol_tag);
    }

    Atom operator-(const Atom& o) const {
        return Atom(pos-o.pos, type);
    }

    double dot(const Atom& o) const
    {
    	return pos.dot(o.pos);
    }

    inline Atom cross(const Atom& B) const
    {
        return Atom(pos.cross(B.pos));
    }

    void rotate(Quat& q)
    {
    	pos.rotate(q);
    }

    inline void rotate(Atom& axis, myFloat angle)
    {
    	pos.rotate(axis.pos, angle);
    }

    void normalise()
    {
        pos.normalise();
    }

    double size() const
    {
    	return pos.size();
    }

    double dist(const Atom& o) const
    {
    	return pos.dist(o.pos);
    }

    double distSQ(const Atom& o) const
    {
    	return pos.distSQ(o.pos);
    }

    bool isNeighbor(const Atom& o, myFloat len = 1.0, myFloat margin = 0.00000001) const
    {
        Atom vec(o.pos.x - pos.x, o.pos.y - pos.y, o.pos.z - pos.z);
        double dist = vec.pos.x*vec.pos.x + vec.pos.y*vec.pos.y + vec.pos.z*vec.pos.z;
        return ( dist < len+margin && dist > len-margin );
    }

    bool isAproxSame(const Atom& o, myFloat approx = 0.000001) const
    {
    	return pos.isAproxSame(o.pos, approx);
    }
};




bool isAproxSame(const myFloat& a, const myFloat& b, myFloat approx = 0.000001) {
    return (a < b+approx && a > b - approx);
}




bool myfunction (Atom i,Atom j) { return (i.pos.size()<j.pos.size()); }




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
    		a *= 1.0/a.pos.dist(zero.pos);
    	}
    }

    void offset(int offs)
    {
        for(Atom& item : (*this))
        {
            item.N += offs;
        }
    }

    //
    //
    // Methods const
    //
    //

    int count_Atoms_of_Type( int atype) const
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




    int count_atoms_of_Mol_tag(int mTag) const
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

    vector<int> get_Mol_Types() const
    {
        vector<int> moltags;
        bool exist = false;
        for(auto& item : (*this))
        {
            exist = false;
            for(auto mtype : moltags)
            {
                if( item.mol_tag == mtype )
                {
                    exist = true;
                }
            }
            if(!exist)
            {
                moltags.push_back(item.mol_tag);
            }
        }
        return moltags;
    }

    int get_Max_Mol_Tag()
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

    Atoms get_molecule(int mol_tag) const
    {
    	Atoms temp;
    	for(auto& a : (*this))
    	{
    		if(mol_tag == a.mol_tag)
    		{
    			temp.push_back(a);
    		}
    	}
    	return temp;
    }

    double get_molecule_radius(int mol_tag) const
    {
    	Atoms molecule = get_molecule(mol_tag);
    	Atom mol_com = center_of_mass(mol_tag);
    	double mol_radius = 0.0;
    	for(auto& m : molecule)
    	{
    		mol_radius += mol_com.dist(m);
    	}
    	return mol_radius / molecule.size();
    }




    bool is_overlap(Atom& b, Force_Field& ff) const
    {
        for(const Atom& a : (*this))
        {
            if( b.dist(a) < ff.get_cutoff(a.type, b.type)*ff.get_cutoff(a.type, b.type) )
            {
                return true;
            }
        }
        return false;
    }

    bool is_overlap(Atoms& other, Force_Field& ff) const
    {
    	//
        // One-All particle overlap
    	//
    	for(Atom& o : other)
    	{
    		if( is_overlap(o, ff) ) // overlap of atom o with this container
    		{
    			return true;
    		}
    	}

    	//
    	// Test if particle is within a hollow molecule
    	//
    	vector<int> molecules = get_Mol_Types();
    	Atom mol_com;
    	double mol_radius;

    	for(auto mol_i : molecules)
    	{
    		mol_com = center_of_mass(mol_i);
    		mol_radius = get_molecule_radius(mol_i);

    		for(Atom& o : other)
    		{
    		 	if( o.dist(mol_com) < mol_radius ) // overlap of atom o with this container
    		 	{
    		    	return true;
    		    }
    		}
    	}

    	return false;
    }




    /**
     * @brief center_of_mass - function computes Center-Of-Mass (COM) of particles with a given mol_tag
     * @param mtag - mol_tag of particles for COM calculation, -1 = all particles regardless of mol_tag
     */
    Atom center_of_mass(int mtag=-1, int start=-1, int stop=-1) const
    {
        int count=0;
        int total=0;
        Atom cm;
        for(const Atom& item : (*this))
        {
            if( (item.mol_tag == mtag || mtag == -1) && total >= start && (total < stop || stop == -1) )
            {
                cm += item;
                ++count;
            }

            if(item.mol_tag == mtag || mtag == -1)
                ++total;
        }
        cm *= 1.0/count;
        return cm;
    }

    //
    //
    // static methods
    //
    //

    static Atom clusterCM(vector<Atom>& cluster)
    {
        Atom cluscm(0.0, 0.0, 0.0);
        for(auto& a : cluster)
        {
            cluscm += a;
        }
        cluscm *= 1.0/cluster.size();

        return cluscm;
    }

    static void clusterRotate(vector<Atom>& cluster, double angle, Tensor_xyz axis) {
        double vc,vs;

        Tensor_xyz cluscm = clusterCM(cluster).pos;

        axis.normalise();

        vc = cos( angle );
        vs = sqrt(1.0 - vc*vc);

        Quat newquat(vc, axis.x*vs, axis.y*vs, axis.z*vs);

        //quatsize=sqrt(newquat.w*newquat.w+newquat.x*newquat.x+newquat.y*newquat.y+newquat.z*newquat.z);

        //shift position to geometrical center
        for(unsigned int i=0; i<cluster.size(); ++i) {
            //shift position to geometrical center
            cluster[i].pos.x -= cluscm.x;
            cluster[i].pos.y -= cluscm.y;
            cluster[i].pos.z -= cluscm.z;
            //do rotation
            cluster[i].pos.rotate(newquat);
            //shift positions back
            cluster[i].pos.x += cluscm.x;
            cluster[i].pos.y += cluscm.y;
            cluster[i].pos.z += cluscm.z;
        }
    }

    static void clusterRotate_random(vector<Atom>& cluster, double max_angle) {
        double vc,vs;
        Tensor_xyz newaxis;

        Tensor_xyz cluscm = clusterCM(cluster).pos;

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
            cluster[i].pos.x -= cluscm.x;
            cluster[i].pos.y -= cluscm.y;
            cluster[i].pos.z -= cluscm.z;
            //do rotation
            cluster[i].pos.rotate(newquat);
            //shift positions back
            cluster[i].pos.x += cluscm.x;
            cluster[i].pos.y += cluscm.y;
            cluster[i].pos.z += cluscm.z;
        }
    }
};

#endif // VECTOR_H
