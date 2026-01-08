#ifndef VECTOR_H
#define VECTOR_H

#include <string>

#include "rng.h"
#include "force_field.h"



class Atom{
public:
    //
    // Lammps atom style full: N molecule-tag atom-type q x y z nx ny nz  (N = # of atoms)
    //
    int N=1;
    int mol_tag = -1;
    int type=0;
    double q=0;
    Tensor_xyz pos = Tensor_xyz(0.0, 0.0, 0.0);
    int nx=0,ny=0,nz=0;
    Tensor_xyz vel = Tensor_xyz(1.001, 1.001, 1.001); // lammps velocities

    //
    // PDB
    //
    // 1-4 ATOM, HETATM, 5-6:empty
    // 7-11: atom serial number
    string atom_serial_N = "    1"; // right, 7 8 9 10 11
    // 13-16:atom name
    string atom_name = " C  "; // left, 13 14 15 16
    // 18-20: Residue name
    string res_name = "THR";  // right
    // 22: chain identifier
    char chain_id = '1';
    // 23-26: Residue sequence number
    int res_seq_N = 1;     // right
    // 27:code for insertion of residues
    char code = ' ';
    // 31-38:x , 39-46:y, 47-54:z :: pos // right
    // 55-60: Occupancy
    double occupancy = 0.0; // right
    // 61-66: Temperature factor
    double temp_factor = 0.0; // right
    // 73-76: Segment identifier (optional)
    string seg_id = "    "; // left
    // 77-78 Element symbol
    string element = " C"; // right
    // 79-80 Charge (optional)
    string charge = "  ";

    //
    // variables for fast loading
    //
    int to_type=-1; // used in dodecahedron to change the type after bonds are generated - design mistake
    double temp_dist = -1.0;
    bool loaded = false; // used for loading lipids

    //
    // Constructors
    //
    Atom() {}
    Atom(Tensor_xyz pos, int type=0): type(type), pos(pos) {}
    Atom(Tensor_xyz pos, int type, int mol_tag): mol_tag(mol_tag), type(type), pos(pos) {}

    Atom(myFloat x, myFloat y, myFloat z, int type=0): type(type), pos(x,y,z) {}
    Atom(myFloat x, myFloat y, myFloat z, myFloat vx, myFloat vy, myFloat vz, int type=0): type(type), pos(x,y,z), vel(vx,vy,vz) {}
    Atom(myFloat x, myFloat y, myFloat z, int type, int mol_tag): mol_tag(mol_tag), type(type), pos(x,y,z) {}

    Atom(int N, Tensor_xyz pos, int type=0): N(N), type(type), pos(pos) {}
    Atom(int N, Tensor_xyz pos, int type, int mol_tag): N(N), mol_tag(mol_tag), type(type), pos(pos) {}

    //
    // Geometric funstions
    //
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

    bool isAproxSame(const Atom& o, myFloat approx = 0.000001) const
    {
        return pos.isAproxSame(o.pos, approx);
    }

    bool isNeighbor(const Atom& o, myFloat len = 1.0, myFloat margin = 0.00000001) const
    {
        Atom vec(o.pos.x - pos.x, o.pos.y - pos.y, o.pos.z - pos.z);
        double dist = vec.pos.x*vec.pos.x + vec.pos.y*vec.pos.y + vec.pos.z*vec.pos.z;
        return ( dist < len+margin && dist > len-margin );
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

    //
    // Operators
    //
    bool operator==(const Atom& o) const
    {
    	return (this->pos == o.pos);
    }

    bool operator!=(const Atom& o) const
    {
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

    //
    // Output
    //
    friend std::ostream& operator<<(std::ostream& os, Atom& a)
    {
        os << a.pos << endl;
        return os;
    }

    //
    // friend functions
    //

};

bool sortN(const Atom& i, const Atom& j) { return i.N < j.N; }
bool isAproxSame(const myFloat& a, const myFloat& b, myFloat approx = 0.000001) { return (a < b+approx && a > b - approx); }
bool myfunction (Atom i,Atom j) { return (i.pos.size()<j.pos.size()); }




class Atoms : public vector< Atom >
{
public:

    void set_frame(vector<Tensor_xyz>& frame)
    {
        this->clear();
        this->resize(frame.size());

        for(unsigned int i=0; i<frame.size(); ++i)
        {
            (*this)[i].pos = frame[i];
        }
    }

    //
    //
    // IO methods
    //
    //
    string toString()
    {
        stringstream ss;
        vector<int> moltags = get_Atom_Types();
        vector<int> types = get_Atom_Types();

        ss << "Beads: " << size() << endl;
        ss << types.size() << " atom types:" << endl;
        for(int atom_type : types)
        {
            ss << "Atom type " << atom_type << " of " << count_Atoms_of_Type(atom_type) << endl;
        }

        ss << moltags.size() << " molTypes:" << endl;
        for(int mol_tag : moltags)
        {
            ss << "Molecule " << mol_tag << " with " << count_atoms_of_Mol_tag(mol_tag) << " atoms" << endl;
        }

        return ss.str();
    }

	//
	//
    // Methods altering items of the container
	//
	//
    /**
     * @brief center - moves structure, so that center_of_mass is (0,0,0)
     * @param mtag
     */
    void center(int mtag=-1)
    {
        if(!empty())
        {
            Atom cm = center_of_mass(mtag);
            cm*=-1.0;
            move( cm );
        }
    }

    void project_to_unit_sphere()
    {
        Atom zero = Atom(0.0,0.0,0.0);

        for(Atom& a : (*this))
        {
            a *= 1.0/a.pos.dist(zero.pos);
        }
    }

	void set_mol_tag(int mtag)
	{
	    for(Atom& item : (*this))
	        item.mol_tag = mtag;
	}

    void add_offset(int offs)
    {
        for(Atom& item : (*this))
        {
            item.N += offs;
        }
    }

    void scale(double scale)
    {
        for(Atom& item : (*this))
            item *= scale;
    }

    void move(Atom move)
    {
        for(Atom& item : (*this))
            item += move;
    }

    //
    // set join with other Atoms container, skip same atoms
    //
    void join(Atoms other)
    {
        bool same = false;
        for(Atom& o : other)
        {
            same = false;
            for(Atom& a : *(this))
            {
                if(a == o)
                    same = true;
            }
            if(!same)
                this->push_back(o);
        }
    }

    //
    //
    // const Methods
    //
    //

    //
    // Getters
    //
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




    //
    // Counters
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





    double min_dist() const
    {
        double dist = 999.9;

        Atom b = this->at(0);
        for(auto& a : *this)
        {
            if(a != b && b.dist(a) < dist)
            {
                dist = b.dist(a);
            }
        }

        return dist;
    }

    bool similar(Atoms& other) const
    {
        for(auto& o : other)
        {
            for(auto& a : *this)
            {
                if(o == a)
                    return true;
            }
        }
        return false;
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

    bool is_overlap(Atom& b, double cutoffSQ=1.0) const
    {
        for(const Atom& a : (*this))
        {
            if( b.distSQ(a) < cutoffSQ )
            {
                return true;
            }
        }
        return false;
    }

    bool is_overlap(Atom& b, double cutoffSQ, vector< vector<int>* >& nei) const
    {
        for(int j=0; j<nei.size(); ++j)
        {
            for(int i : (*nei[j]))
            {
                if( b.distSQ( (*this)[i] ) < cutoffSQ )
                {
                    return true;
                }
            }
        }
        return false;
    }

    bool is_within_hollow(Atoms& other) const
    {
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

    bool is_overlap(Atoms& other, double cutoff=1.0, bool check_hollow=false) const
    {
        for(Atom& o : other)
        {
            if( is_overlap(o, cutoff*cutoff) ) // overlap of atom o with this container
            {
                return true;
            }
        }

        if(check_hollow)
            return is_within_hollow(other);

        return false;
    }

    bool is_overlap(Atoms& other, double cutoff, vector< vector<int>* >& nei) const
    {
        for(Atom& o : other)
        {
            if( is_overlap(o, cutoff*cutoff, nei) ) // overlap of atom o with this container
            {
                return true;
            }
        }

        return false;
    }

    bool is_overlap(Atoms& other, Force_Field& ff, bool check_hollow=false) const
    {
    	for(Atom& o : other)
    	{
    		if( is_overlap(o, ff) ) // overlap of atom o with this container
    		{
    			return true;
    		}
    	}

        if(check_hollow)
            return is_within_hollow(other);

        return false;
    }




    Atom get_center_of_mass() const
    {
        return this->center_of_mass(-1,-1,-1);
    }

    Atom center_of_mass(int mtag) const
    {
        return this->center_of_mass(mtag, -1,-1);
    }

    Atom center_of_mass(int start, int stop) const
    {
        return this->center_of_mass(-1, start, stop);
    }

    /**
     * @brief center_of_mass - function computes Center-Of-Mass (COM) of particles with a given mol_tag
     * @param mtag - mol_tag of particles for COM calculation, -1 = all particles regardless of mol_tag
     */
    Atom center_of_mass(int mtag, int start, int stop) const
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

    Atoms within(Atom& start, double radius)
    {
        Atoms select;

        for(auto& a : *(this))
        {
            if( start.dist(a) < radius )
            {
                select.push_back(a);
            }
        }

        return select;
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
