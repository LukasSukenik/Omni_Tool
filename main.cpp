#include "rng.h"
#include "icosahedron.h"
#include "oblatespheroid.h"
#include "sphere.h"
#include "tennisball.h"
#include "spherepatch.h"
#include "pentamer.h"
#include "dodecahedron.h"
#include "chain.h"
#include "slab.h"
#include "globular_sphere.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>

using namespace std;




/**
 * @brief Structure_Container
 * map classes, each generating a different structure
 */
class Structure_Container : public map<string, Particles*> // @suppress("Invalid template argument")
{
public:
    Structure_Container()
    {
        //
        // Add structures
        //
        (*this)[Empty_Particle::keyword] = new Empty_Particle();
        (*this)[Icosahedron<Surface>::keyword] = new Icosahedron<Surface>();
        (*this)[Sphere::keyword] = new Sphere();
        (*this)[TennisBall::keyword] = new TennisBall();
        (*this)[OblateSpheroid::keyword] = new OblateSpheroid();

        (*this)[SpherePatch::keyword] = new SpherePatch();
        (*this)[Pentamer<PentaSurface>::keyword] = new Pentamer<PentaSurface>();
        (*this)[Dodecahedron::keyword] = new Dodecahedron();
        (*this)[Chain::keyword] = new Chain();
        (*this)[Slab::keyword] = new Slab();

        (*this)[NettedSlab::keyword] = new NettedSlab();
        (*this)[SphereJanus::keyword] = new SphereJanus();
        (*this)[Cow::keyword] = new Cow();
        (*this)[Globular_Sphere::keyword] = new Globular_Sphere();
        (*this)[Monomer::keyword] = new Monomer();
    }

    ~Structure_Container()
    {
        for (auto const& [key, val] : (*this)) // @suppress("Symbol is not resolved")
        {
             delete val;
        }
    }

	void helpMessage()
	{
		cout << "Specify input files, run $ ./ico filename_1 filename_2 ...\n\n" << endl;
		cout << "Example input files:\n" << endl;

		for (auto& [key, val] : *this) // @suppress("Symbol is not resolved")
		{
		    std::cout << (*val).help() << std::endl; // @suppress("Method cannot be resolved") // @suppress("Invalid overload")
		}
	}
};

void do_analysis()
{
	Atom set_a[3];
	set_a[0] = Atom(0,0,0);
	set_a[1] = Atom(0,-1.2508,0);
	set_a[2] = Atom(0, 1.2508,0);

	Atom b = Atom(0,0,0);

	Force_Field ff;
	ff.lj[0] = LJ(0, 1.0, 1.75652, 1.97162);
	ff.cos[0] = CosSQ(0, 1.75652, 1.97162, 0.5);

	int size=600;

	double step=0.01;
	double e,f;
	Tensor_xyz f_vec, f_sum;
	bool inside = false;
	bool inside2 = false;
	double min=0,max=0,value=0;

	for(int i=-size; i<=size; ++i)
	{
		for(int j=-size; j<=size; ++j)
		{
			b.pos = Tensor_xyz(step*i, step*j, 0.0);
			e=0.0;
			f_sum = Tensor_xyz();
			inside = false;
			inside2 = false;

			for(auto a : set_a)
			{
				e += ff.energy(a.dist(b), a.type, b.type);
				f = ff.force(a.dist(b), a.type, b.type);
				f_vec = b.pos-a.pos;
				f_vec.normalise();
				f_vec *= f;
				f_sum += f_vec;

				if(a.dist(b) < 1.97162) inside = true;
			    if(a.dist(b) < 0.5*1.97162) inside2 = true;
			}
			value=abs(f_sum.y); // f_sum.size() e
			cout << ( (inside) ? ( (inside2) ? -0.02 : -0.01 ) : value ) << " ";

			if(!inside) // calc min, max
			{
				if(value < min) min = value;
				if(value > max) max = value;
			}
		}
		cout << endl;
	}

	cerr << min << " " << max << endl;

	/*std::ofstream outFile("vars.sh");
	outFile << "bead=" << min-0.02 << endl;
	outFile << "sigma=" << min-0.01 << endl;
	outFile << "min=" << max*0.5 << endl;
	outFile << "mmm=" << max*0.9 << endl;
	outFile << "max=" << max << endl;
	outFile.close();*/
}


int main(int argc, char* argv[]) // // $num of beads per edge, box dimensions X(same as Y) $beg $end, $position in Z, $offset
{
    Data data;
    Structure_Container structure;

    //
    // Input safeguard
    //
    if( argc == 1 || strcmp(argv[1], "-h") == 0 )
    {
        structure.helpMessage();
        exit(1);
    }

    //
    // Loop over command line arguments: input file names
    //
    for(int i=1; i<argc; ++i)
    {
        data.loadInput(argv[i]); // Load the input files specified as command line arguments
        cerr << data.in.toString() << endl; // Report what was loaded

        //
        // Load a particle if specified
        //
        if( data.isDefined() )
        {
            data.load(data.in.infile);    // Load Data from infile
            data.scale(data.in.scale);    // Rescale atom positions
            data.move(data.in.com_pos);   // Move entire system by vector

            if(data.in.is_mol_tag())
            	data.set_mol_tag(data.in.mol_tag); // Change mol_tag of all particles to one set by input
            if(data.in.is_mtag_12())
                data.align(data.in.mtag_1, data.in.mtag_2); // align mol_tag particles in z axis and XY plane

            // if generating into an existing structure that you did not load, give the number of particles as offset
            data.offset(data.all_beads.size());

            if( data.in.fit )
                data.fit();
            if( data.in.center )
                data.center();

            data.add(); // temp_data to all_data
        }

        //
        // generate a particle when no input particle is specified
        //
        if( ! data.isDefined() )
        {
            if(structure.count(data.in.gen_structure) > 0) // @suppress("Method cannot be resolved")
            {
                cerr << "Generating: " << structure[ data.in.gen_structure ]->name << endl; // @suppress("Field cannot be resolved") // @suppress("Invalid overload")
				structure[ data.in.gen_structure ]->generate( data ); // @suppress("Method cannot be resolved")
                structure[ data.in.gen_structure ]->scale( data.in.scale ); // @suppress("Method cannot be resolved")
                structure[ data.in.gen_structure ]->move( data.in.com_pos ); // @suppress("Method cannot be resolved")
                structure[ data.in.gen_structure ]->populate( data ); // @suppress("Method cannot be resolved")
                structure[ data.in.gen_structure ]->add(data); // particle data // @suppress("Method cannot be resolved")
            }
            else
            {
            	cerr << "not found keyword " << data.in.gen_structure << endl;
            	exit(2);
            }
        }
    }

    //
    // Analyze stuff
    //
    if(false)
    	do_analysis();

    vector<string> coeff;
    vector<double> dist = data.createBondGroups(coeff);
    //data.roundBondDist(dist);
    //data.removeDuplicateBond();

    data.print();

    //
    // Report on the Final structure
    //
    if(data.all_beads.empty())
    {
        cerr << "!!! Nothing generated !!!" << endl;
    }
    else
    {
        cerr << data.toString() << endl;
    }

    return 0;
}
