#ifndef SYSTEM_CONTAINER_H
#define SYSTEM_CONTAINER_H

#include <iostream>
#include <map>

#include "system_base.h"
#include "virus_pseudot3.h"
#include "membrane.h"
#include "lipid_nanoparticle.h"
#include "vesicle.h"

using namespace std;

/**
 * @brief Structure_Container
 * map classes, each generating a different structure
 */
class System_Container : public map<string, System_Base*>
{
public:
    System_Container()
    {
        //
        // Add systems
        //
        (*this)[Virus_pseudoT3::keyword] = new Virus_pseudoT3();
        (*this)[Vesicle::keyword] = new Vesicle();
        (*this)[Lipid_Nanoparticle::keyword] = new Lipid_Nanoparticle();
        //(*this)[Icosahedron<Surface>::keyword] = new Icosahedron<Surface>();
    }

    ~System_Container()
    {
        for (auto const& [key, val] : (*this))
        {
             delete val;
        }
    }

    void helpMessage()
    {
        cout << "Specify input files, run $ ./omni_tool filename_1 filename_2 ...\n\n" << endl;
        cout << "Example input files:\n" << endl;

        for (auto& [key, val] : *this)
        {
            cout << (*val).help() << endl;
        }
        exit(1);
    }
};

#endif // SYSTEM_CONTAINER_H
