#ifndef MEHRNOOSH_SPHERE_H
#define MEHRNOOSH_SPHERE_H

#include "sphere.h"

class Globular_Sphere : public Sphere
{
public:
    inline static const string keyword = "globular_sphere";
    const string name = "globular_sphere";

    Globular_Sphere() : Sphere("Globular_Sphere") {}

    void generate( Data& data )
    {
        Atoms base_sphere;
        Atoms patch_sphere;
        Atoms full_sphere;

    	// generate low count pseudoatom sphere, for repuls only interation
        if(data.in.subdiv_beads != -1 && data.in.subdiv_lig != -1)
        {
        	base_sphere = icosphere(data.in.subdiv_beads, data.in.ff.lj[1].type, data.in.mol_tag);
      	}
       	if(data.in.num_of_beads != -1 && data.in.num_lig != -1)
       	{
       		fibonacci_sphere( base_sphere, data.in.num_of_beads, data.in.ff.lj[1].type, data.in.mol_tag);
       	}


    	// generate sphere for patches
    	if(data.in.subdiv_beads != -1 && data.in.subdiv_lig != -1)
    	{
    		patch_sphere = icosphere(data.in.subdiv_lig, 0, data.in.mol_tag);
    	}
    	if(data.in.num_of_beads != -1 && data.in.num_lig != -1)
    	{
    		fibonacci_sphere( patch_sphere, data.in.num_lig, -1, data.in.mol_tag);
    	}

    	// scale patch_sphere by 1.0 + 1st sphere pseudoatom sigma - patch pseudoatom sigma
    	patch_sphere.scale(1.0 + (data.in.ff.lj[1].sigma - data.in.ff.lj[2].sigma) / data.in.scale);

    	// create the final particle based on defined patches
    	for(auto& atom : patch_sphere)
    	{
    		for(auto& patch : data.in.patches)
    	    {
    			if( (patch.x*(atom.x - patch.vx)) + (patch.y*(atom.y - patch.vy)) + (patch.z*(atom.z - patch.vz)) > 0 + (data.in.ff.lj[1].sigma - data.in.ff.lj[2].sigma) / data.in.scale )
    	        {
    				atom.type=patch.type;
    	        	full_sphere.push_back(atom);
    	        }
    	    }
    	}

    	bool add = true;
    	for(auto& atom : base_sphere)
    	{
    		add = true;
    		for(auto& patch : data.in.patches)
    		{
    	     	if( (patch.x*(atom.x - patch.vx)) + (patch.y*(atom.y - patch.vy)) + (patch.z*(atom.z - patch.vz)) > 0 )
    	       	{
    	        	add = false;
    	     	}
    	  	}
    	    if(add)
    	    {
    	    	full_sphere.push_back(atom);
    	    }
    	}

        if( data.in.subdiv_beads == -1 && data.in.subdiv_lig == -1 &&
			data.in.num_of_beads == -1 && data.in.num_lig == -1)
        {
        	cerr << "Specify either:" << endl;
        	cerr << "Number_of_beads: number_of_beads_of_a_sphere" << endl;
        	cerr << "Number_of_ligands: number_of_beads_of_a_sphere" << endl;
        	cerr << "or:" << endl;
        	cerr << "Subdiv_of_beads: 1-5" << endl;
        	cerr << "Subdiv_of_ligands: 1-5" << endl;
        }

        // insert the final structure to beads
        beads.insert(beads.end(), full_sphere.begin(), full_sphere.end());
    }

protected:

    Atoms icosahedron(int type, int mol_tag)
    {
    	Atoms ico;

    	double a = 1.0;
    	double b = 1.0 / sqrt(5.0);
    	double c = 2.0 / sqrt(5.0);
    	double d = (5.0-sqrt(5.0))/10.0;
    	double e = (-5.0-sqrt(5.0))/10.0;
    	double f = sqrt(-e);
    	double g = sqrt(d);

    	// add vertices
    	ico.push_back( Atom( 0.0, 0.0, -a, type, mol_tag) );
    	ico.push_back( Atom( 0.0, 0.0,  a, type, mol_tag) );
    	ico.push_back( Atom(   c, 0.0,  b, type, mol_tag) );
    	ico.push_back( Atom(   d,   f,  b, type, mol_tag) );
    	ico.push_back( Atom(   e,   g,  b, type, mol_tag) );
    	ico.push_back( Atom(   e,  -g,  b, type, mol_tag) );
    	ico.push_back( Atom(   d,  -f,  b, type, mol_tag) );
    	ico.push_back( Atom(  -c, 0.0, -b, type, mol_tag) );
    	ico.push_back( Atom(  -d,  -f, -b, type, mol_tag) );
    	ico.push_back( Atom(  -e,  -g, -b, type, mol_tag) );
    	ico.push_back( Atom(  -e,   g, -b, type, mol_tag) );
    	ico.push_back( Atom(  -d,   f, -b, type, mol_tag) );

    	return ico;
    }

    void subdivide(Atoms& ico)
    {
    	double min_dist = 999.9;
    	Atoms add;
    	Atom a,b,c;

    	for(int i=0; i<ico.size(); ++i)
    	{
    		for(int j=0; j<i; ++j)
    		{
    			a = ico[i];
    			b = ico[j];
    			if(a != b && a.dist(b) < min_dist)
    			{
    				min_dist = a.dist(b);
    			}
    		}
    	}

    	for(int i=0; i<ico.size(); ++i)
    	{
    	 	for(int j=0; j<i; ++j)
    	   	{
    	 		a = ico[i];
    	 		b = ico[j];
    	    	if(a != b && a.dist(b) > min_dist-0.01 && a.dist(b) < min_dist+0.01 )
    	    	{
    	    		c = (a + b) * 0.5;
    	    		add.push_back(c);
    	    	}
    	    }
    	}

    	ico.insert( ico.end(), add.begin(), add.end() );
    }

    Atoms icosphere(int subdiv, int type, int mol_tag)
    {
    	Atoms ico = icosahedron(type, mol_tag);

    	for (size_t i = 1; i < subdiv; ++i)
    	{
    	    subdivide(ico);
    	}
    	ico.project_to_unit_sphere();
        return ico;
    }
};

#endif // MEHRNOOSH_SPHERE_H
