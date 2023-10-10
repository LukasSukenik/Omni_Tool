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

using namespace std;

void helpMessage(vector<Particles*>& nano);

int main(int argc, char* argv[]) // // $num of beads per edge, box dimensions X(same as Y) $beg $end, $position in Z, $offset
{
    Data data;
    vector<Particles*> nano;
    nano.push_back( new Empty_Particle() );                      // 1
    nano.push_back( new Icosahedron<Surface>("Icosahedron") );                      // 2
    nano.push_back( new Sphere() );                                                 // 3
    nano.push_back( new TennisBall() );                                             // 4
    nano.push_back( new OblateSpheroid() );                                         // 5
    nano.push_back( new SpherePatch() );                                            // 6
    nano.push_back( new Pentamer<PentaSurface>("Icosahedron + PentaSurface" ) );    // 7
    nano.push_back( new Dodecahedron("Dodecahedron") );                             // 8
    nano.push_back( new Chain() );                                                  // 9
    nano.push_back( new Slab() );                                                   // 10
    nano.push_back( new NettedSlab());                                              // 11
    nano.push_back( new SphereJanus());                                             // 12
    nano.push_back( new Cow());                                                     // 13

    if( argc == 1 || strcmp(argv[1], "-h") == 0 ) {
        cout << "No input file specified" << endl;
        helpMessage(nano);
        exit(1);
    }

    for(int i=1; i<argc; ++i)
    {
        data.loadInput(argv[i]);

        cerr << data.in.toString() << endl;

        if( data.isDefined() )
        {
            data.load(data.in.infile);      // Load Data from file "data.in.infile"
            if(data.in.isScale())
                data.rescale(data.in.scale);    // Rescale atom positions by data.in.scale - usually 1
            if(data.in.isCOM_pos())
                data.move(data.in.com_pos);     // Move by vector defined in input file
            if(data.in.is_mol_tag())
                data.mol_tag(data.in.mol_tag);  // Change mol_tag of all particles to one set by input
            if(data.in.is_mtag_12())
                data.align(data.in.mtag_1, data.in.mtag_2); // align mol_tag particles in z axis and XY plane

            data.offset(data.all_beads.size());

            if( data.in.fit )
                data.fit();
            if( data.in.center )
                data.center();

            data.add(); // temp_data to all_data
        }
        else
        {
            if(data.in.nano_type < 1 || data.in.nano_type > nano.size()) {
                cerr << "Bad nanoparticle type selected" << endl;
                exit(2);
            }

            cerr << "Generating: " << nano[data.in.nano_type-1]->name << endl;
            nano[data.in.nano_type-1]->generate( data );
            nano[data.in.nano_type-1]->rescale( data.in.scale );
            nano[data.in.nano_type-1]->move( data.in.com_pos );
            nano[data.in.nano_type-1]->add(data); // particle data
        }
    }

    vector<string> coeff;
    vector<double> dist = data.createBondGroups(coeff);
    //data.roundBondDist(dist);
    //data.removeDuplicateBond();

    if( data.in.out_type == 1)
    {
        data.printLammps();
        cerr << "Lammps out" << endl;
    }
    if( data.in.out_type == 0)
    {
        data.printXYZ();
        cerr << "XYZ out" << endl;
    }
    if( data.in.out_type == 2)
    {
    	data.printPDB();
    	cerr << "PDB out" << endl;
    }

    for(int i=0; i<nano.size(); ++i)
        delete nano[i];

    cerr << data.toString() << endl;

    return 0;
}




void helpMessage(vector<Particles*>& nano)
{
    cout << "./ico filename_1 filename_2" << endl;

    cout << "\nGenerated structure types keyword is \"Particle_type: number\"" << endl;

    for(int i=0; i< nano.size(); ++i)
    {
    	cout << i+1 << " is " << nano[i]->name << endl;
    }

    cout << "Keywords:" << endl;
    cout << "Ouput/Input category:" << endl;
    cout << "Output_type: 1 = Lammps" << endl;
    cout << "             0 = XYZ" << endl;
    cout << "Lammps_offset: integer" << endl;
    cout << " - offset the generated structure for manual insertion into another lammps structure file" << endl;
    cout << "Load_file: filename" << endl;


    cout << "Num_of_beads: integer" << endl;
    cout << " - Particle_type: 1,6,7 = number of beads per edge" << endl;
    cout << " - other Particle_type = number of beads for entire nanoparticle/structure" << endl;
    cout << "Scale: floating_point_number - nanoparticle radius" << endl;
    cout << "c: float " << endl;
    cout << " - Particle_type: 4 = oblate spheroid < 1, prolate spheroid > 1, 1.0 - ERROR, not defined" << endl;
    cout << " - Particle_type: 3 = width of patch (0.0 to 2.0)" << endl;
    cout << "Number_of_ligands: integer" << endl;
    cout << "Mol_tag: integer" << endl;
    cout << " - change mol_tag of generated/loaded structure " << endl;
    cout << "Atom_type: integer" << endl;
    cout << " - Atom_type of generated structure, if structure has more atom_types they are incremented from provided value " << endl;
    cout << "Janus: float float float" << endl;

    cout << "\nPosition/Box properties category:" << endl;
    cout << "Box: float float float float float float" << endl;
    cout << "position_shift: float float float" << endl;
    cout << "Center = centers particles at 0.0" << endl;
    cout << "Align: mol_tag_1 mol_tag_2 integer" << endl;
    cout << " - align mol_tag_1 with x-axis" << endl;
    cout << " - center mol_tag_2 around z axis" << endl;
    cout << "Impact_vector: float float float" << endl;
    cout << "Fit" << endl;
    cout << " - positions the loaded/generated structure next to previosly generated/loadedd structure" << endl;
    cout << " - used for ideal collision position of two liposomes and a nanoparticle" << endl;

    cout << "\nForce-Field category:" << endl;
    cout << "Beads_lj/cut:" << endl;

    cout << "\nSeed: intege = radom generator" << endl;
}
