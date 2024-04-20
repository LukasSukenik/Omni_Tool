#ifndef SPHERE_H
#define SPHERE_H

#include "particle.h"

using namespace std;

class Sphere: public Particles {
public:
    inline static const string keyword = "sphere";

    const string name = "sphere";

    Sphere() : Particles("sphere") {}
    Sphere(string str) : Particles(str) {}

    void generate( Data& data )
    {
        vector<Atom> ligand;
        int nano_start = beads.size();
        data.in.c=1.0;

        fibonacci_sphere( beads, data.in.num_of_beads, typeNano);
        fibonacci_sphere_z_distrib_linear( ligand, data.in.num_lig, data.in.c, typeTemp); // second fib.

        for(auto& patch : data.in.patches)
        {
            gen_ligands_2( data, ligand, patch, typeNano); // find closest spheres on first fib. sphere and change their type
        }

        beads.erase(beads.begin()+data.in.num_of_beads, beads.end()); // erase second fib sphere

        int nano_end = beads.size();
        for(auto& atom : ligand)
        {
        	atom.mol_tag = data.in.mol_tag;
        }
    }


    void fibonacci_sphere(vector<Atom>& container, int samples, int type, int mol_tag=0)
    {
        const double PI = 3.141592653589793;
        double offset = 2.0/samples;
        double increment = PI * (3.0 - sqrt(5.0));
        double x,y,z,r,phi;

        for(int i=0; i<samples; ++i) {
            z = ((i * offset) - 1) + (offset / 2);
            r = sqrt(1 - z*z);

            phi = i * increment;

            x = cos(phi) * r;
            y = sin(phi) * r;

            container.push_back(Atom(x,y,z,type, mol_tag));
        }
    }

    void move(Atom move)
    {
        for(Atom& item : this->beads)
            item += move;
    }

    void rescale(double rescale)
    {
        for(Atom& item : this->beads)
            item *= rescale;
    }

    string help()
    {
    	stringstream ss;

    	ss << "Particle_type: sphere\n";
    	ss << "Output_type: pdb # other keywords: xyz pdb lammps_full\n";
		ss << "Num_of_beads: 500\n";
		ss << "Scale: 5.0\n";
		ss << "Number_of_ligands: 50\n";
        ss << "Mol_tag: 2\n";

    	return ss.str();
    }

protected:    
    void gen_ligands( Data& data, vector<Atom>& ligand, Atom patch, int type_from, int type_temp=-1)
    {
        for(auto& lig : ligand)
        {
            if(lig.pos.x < patch.pos.x && lig.pos.x > patch.vel.x &&
               lig.pos.y < patch.pos.y*data.in.c && lig.pos.y > patch.vel.y*data.in.c &&
               lig.pos.z < patch.pos.z && lig.pos.z > patch.vel.z)
            {
                Atom* select = &beads[0]; // Stupid C++, for some reason reference dont work
                for(auto& item : beads)
                {
                    if(select->dist(lig) > item.dist(lig) && item.type == type_from )
                    {
                        select = &item;
                    }
                }
                if( select->type == type_from )
                    select->type = patch.type;
            }
        }
    }

    void gen_ligands_2( Data& data, vector<Atom>& ligand, Atom patch, int type_from, int type_temp=-1)
    {
        for(auto& lig : ligand)
        {
            // plane equation a*x-x_coord + ... > 0
            if( (patch.pos.x*(lig.pos.x - patch.vel.x)) + (patch.pos.y*(lig.pos.y - patch.vel.y)) + (patch.pos.z*(lig.pos.z - patch.vel.z)) > 0 )
            {
                Atom* select = &beads[0];
                for(auto& item : beads)
                {
                    if(select->dist(lig) > item.dist(lig) && item.type == type_from )
                    {
                        select = &item;
                    }
                }
                if( select->type == type_from )
                    select->type = patch.type;
            }
        }
    }

    /**
     * @brief fibonacci_sphere_z_distrib_linear - We are transforming z axis coordinates where we generate particles. In sphere its linear - slope 1
     *                                          - We want x from -1 to 0
     * @param samples
     * @param type
     */
    void fibonacci_sphere_z_distrib_linear(vector<Atom>& container, int samples, double S, int type)
    {
        const double PI = 3.141592653589793;
        double offset = 2.0/samples;
        double increment = PI * (3.0 - sqrt(5.0));
        double x,y,z,r,phi;

        for(int i=0; i<samples; ++i) {

            z = ((i * offset) - 1) + (offset / 2); // z linearly distributed from (-1, -1) to (-1, -1)
            z = transform(z, S);
            r = sqrt(1 - z*z);

            phi = i * increment;

            x = cos(phi) * r;
            y = sin(phi) * r;

            container.push_back(Atom(x,y,z,type));
        }
    }


    void testDistribution()
    {
        int size = 20;
        int hist[size] = {0};
        for(int i=0; i<beads.size(); ++i) {
            if(beads[i].type == typeLig) {
                cout << beads[i].pos.z << endl;
                ++hist[ (int)(beads[i].pos.z*20) ];
            }
        }
        for(int i=0; i< size; ++i) {
            //cout << hist[i] << endl;
        }
        exit(1);
    }

private:
    double transform(double z, double param)
    {
        if(z < 0)
            return -( param * pow(fabs(z), 0.5)) - (1-param);
        else
            return ( param * pow(fabs(z), 0.5)) + (1-param);
    }
};

class SphereJanus : public Sphere
{
public:
    inline static const string keyword = "sphere_janus";
    const string name = "sphere_janus";

    SphereJanus() : Sphere("sphere_janus") {}

    void generate( Data& data ) {

        // generate nano
        int nano_start = beads.size();
        fibonacci_sphere(beads, data.in.num_of_beads, typeNano);
        int nano_end = beads.size();

        // change type nano to type lig based on placement
        for(int i=nano_start; i<nano_end; ++i) {
            if(beads[i].pos.x > 0.0) {
                beads[i].type = typeLig;
            }
        }

        for(int i=nano_start; i<nano_end; ++i) {
            beads[i] = ((beads[i]*data.in.scale) + data.in.com_pos);
            beads[i].mol_tag = 2;
        }

        /*nano_start = beads.size();
        fibonacci_sphere(data.num_of_beads, typeNano);
        nano_end = beads.size();

        // change type nano to type lig based on placement
        for(int i=nano_start; i<nano_end; ++i) {
            if(beads[i].x < 0.0) {
                beads[i].type = typeLig;
            }
        }

        data.com_pos.x = 50.0-data.com_pos.x;

        for(int i=nano_start; i<nano_end; ++i) {
            beads[i] = ((beads[i]*data.scale) + data.com_pos);
            beads[i].mol_tag = 3;
        }*/
    }
};

#endif // SPHERE_H
