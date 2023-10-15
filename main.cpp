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
class Structure_Container : public map<string, Particles*>
{
public:
	void helpMessage()
	{
		cout << "Specify input files, run $ ./ico filename_1 filename_2 ...\n\n" << endl;
		cout << "Example input files:\n" << endl;

		for (auto& [key, val] : *this)
		{
		    std::cout << (*val).help() << std::endl;
		}
	}
};




int main(int argc, char* argv[]) // // $num of beads per edge, box dimensions X(same as Y) $beg $end, $position in Z, $offset
{
    Data data;

    Structure_Container structure;

    //
    // Adding a structure generating class
    //
    structure[Empty_Particle::keyword] = new Empty_Particle();
    structure[Icosahedron<Surface>::keyword] = new Icosahedron<Surface>();
    structure[Sphere::keyword] = new Sphere();
    structure[TennisBall::keyword] = new TennisBall();
    structure[OblateSpheroid::keyword] = new OblateSpheroid();

    structure[SpherePatch::keyword] = new SpherePatch();
    structure[Pentamer<PentaSurface>::keyword] = new Pentamer<PentaSurface>();
    structure[Dodecahedron::keyword] = new Dodecahedron();
    structure[Chain::keyword] = new Chain();
    structure[Slab::keyword] = new Slab();

    structure[NettedSlab::keyword] = new NettedSlab();
    structure[SphereJanus::keyword] = new SphereJanus();
    structure[Cow::keyword] = new Cow();

    //
    // Input safeguard
    //
    if( argc == 1 || strcmp(argv[1], "-h") == 0 ) {
        structure.helpMessage();
        exit(1);
    }

    //
    // Loop over command line arguments
    //
    for(int i=1; i<argc; ++i)
    {
        //
        // Load the input file [i]
        //
        data.loadInput(argv[i]);

        //
        // Report what was loaded
        //
        cerr << data.in.toString() << endl;

        //
        // Load a particle if specified
        //
        if( data.isDefined() )
        {
            data.load(data.in.infile);      // Load Data from file "data.in.infile"
            data.scale(data.in.scale);    // Rescale atom positions by data.in.scale
            data.move(data.in.com_pos);     // Move by vector defined in input file

            if(data.in.is_mol_tag())
                data.mol_tag(data.in.mol_tag);  // Change mol_tag of all particles to one set by input
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
            if(structure.count(data.in.gen_structure) > 0)
            {
                cerr << "Generating: " << structure[ data.in.gen_structure ]->name << endl;
				structure[ data.in.gen_structure ]->generate( data );
                structure[ data.in.gen_structure ]->scale( data.in.scale );
                structure[ data.in.gen_structure ]->move( data.in.com_pos );
                structure[ data.in.gen_structure ]->populate( data );
                structure[ data.in.gen_structure ]->add(data); // particle data
            }
            else
            {
            	cerr << "not found keyword " << data.in.gen_structure << endl;
            	exit(2);
            }
        }
    }

    vector<string> coeff;
    vector<double> dist = data.createBondGroups(coeff);
    //data.roundBondDist(dist);
    //data.removeDuplicateBond();

    data.print();

    //
    // free memory
    //
    for (auto const& [key, val] : structure)
    {
         delete val;
    }

    //
    // Report on the Final structure
    //
    cerr << data.toString() << endl;

    return 0;
}
