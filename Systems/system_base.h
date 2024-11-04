#ifndef SYSTEM_BASE_H
#define SYSTEM_BASE_H

#include <sstream>
#include <string>

#include "data.h"

using namespace std;

/**
 * @brief The System_Base class - A virtual class serving as an interface for system generation
 */
class System_Base
{
public:
    inline static const string keyword = "Abstract class of System";
    const string name = "Abstract class of System";

    System_Base() {}
    System_Base(string name) : name(name) {}

    /**
     * @brief generate - Generate particle
     * - Particle need to be generated in scale and in position (data.in.scale and data.in.com_pos)
     * - Position overlap check to already existing particles in data.all_beads
     *
     * @param data
     */
    virtual void generate( Data& data )=0;

    virtual string help()
    {
        stringstream ss;
        ss << "Abstract class System\n";
        ss << "Contains functions intended for inheritance\n";
        ss << "Does not generate anything\n";
        return ss.str();
    }
};

#endif // SYSTEM_BASE_H
