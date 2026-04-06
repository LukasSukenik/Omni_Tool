#ifndef COW_H
#define COW_H

#include <fstream>
#include <sstream>
#include "atom.h"
#include "particle.h"

using namespace std;

class Cow : public Particle
{
public:
    inline static const string keyword = "cow";
    const string name = "cow";

    Cow() : Particle("cow")  {}

    string help()
    {
        stringstream ss;
        ss << "Particle_type: cow\n";
        return ss.str();
    }


    void generate( Data& data  )
    {
        int type=2;
        int mol_tag=2;
        int size=0;
        char line[256];
        string temp;
        double x,y,z;
        std::ifstream fs( "cow3.xyz" );

        if (!fs.is_open())
        {
            std::cout << "failed to open " << "cow3.xyz" << '\n';
        }
        else
        {
            fs >> size;
            fs.getline(line, 256);
            for(int i=0; i<size; ++i) {
                fs >> temp >> x >> y >> z;
                beads.push_back( Atom(x, y, z, type, mol_tag) );
            }
        }
        fs.close();
    }
};

#endif // COW_H
