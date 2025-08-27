#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <map>

#include "particle_container.h"
#include "system_container.h"

using namespace std;

void do_analysis();




class Version
{
public:
    int v=3;

    Version()
    {
        splash();
    }

    void splash()
    {
        cerr << "*******************" << endl;
        cerr << "*                 *" << endl;
        cerr << "* Omni Tool :: v" << v << " *" << endl;
        cerr << "*                 *" << endl;
        cerr << "*******************" << endl;
    }
};




int main(int argc, char* argv[])
{
    Version v;

    Data data;
    Particle_Container particles;
    System_Container systems;

    //
    // Input safeguard
    //
    if( argc == 1 || strcmp(argv[1], "-h") == 0 )
    {
        particles.helpMessage();
    }

    //
    // Loop over command line arguments: input file names
    //
    for(int i=1; i<argc; ++i)
    {
        data.load_input(argv[i]); // Load the input files specified as command line arguments
        cerr << data.in.toString() << endl; // Report what was loaded

        //
        // Load a particle if specified
        //
        if( data.is_load_file() )
        {
            data.load_data(data.in.file_structure); // Load Data from infile
            data.modify();                  // Modify the loaded Data
            // analyze
        }

        //
        // generate a particle if specified
        //
        if( data.is_particle_gen() )
        {
            if(particles.count(data.in.gen_structure) > 0)
            {
                cerr << "Generating particle: " << particles[ data.in.gen_structure ]->name << endl;
                particles[ data.in.gen_structure ]->generate( data );
                particles[ data.in.gen_structure ]->modify( data );
                particles[ data.in.gen_structure ]->populate( data );
                particles[ data.in.gen_structure ]->make_persistent(data); // particle data
            }
            else
            {
                cerr << "main.cpp particle keyword not found " << data.in.gen_structure << endl;
            	exit(2);
            }
        }

        //
        // execute a system-wide method if specified
        //
        if( data.is_system() )
        {
            if(systems.count(data.in.system_type) > 0)
            {
                cerr << "Loading system: " << systems[ data.in.system_type ]->name << endl;
                systems[ data.in.system_type ]->execute( data );
            }
            else
            {
                cerr << "main.cpp system :keyword: not found " << data.in.system_type << endl;
                exit(2);
            }
        }
    }

    data.print();

    return 0;
}




void do_analysis() // capsid stabilit
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
