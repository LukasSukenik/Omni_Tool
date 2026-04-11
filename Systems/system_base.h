#ifndef SYSTEM_BASE_H
#define SYSTEM_BASE_H

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

    const double PI = 3.141592653589793;
    const double deg_to_rad = 0.0174532925199;
    const double rad_to_deg = 57.2957795131;

    System_Base() {}
    System_Base(string name) : name(name) {}

    /**
     * @brief generate - Generate particle
     * - Particle need to be generated in scale and in position (data.in.scale and data.in.com_pos)
     * - Position overlap check to already existing particles in data.all_beads
     *
     * @param data
     */
    virtual void execute( Data& data )=0;

    virtual string help()
    {
        return "";
    }
};

#endif // SYSTEM_BASE_H
