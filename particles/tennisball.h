#ifndef TENNISBALL_H
#define TENNISBALL_H

#include "sphere.h"

using namespace std;

class TennisBall: public Sphere {
private:
    int angle = 10;

public:
    inline static const string keyword = "tennisball";
    const string name = "tennisball";

    TennisBall() : Sphere("tennisball") {}

    void generate( Data& data )
    {
        vector<Atom> ligand;
        fibonacci_sphere(beads, data.in.num_of_beads, typeNano);
        fibonacci_sphere_z_distrib_linear(ligand, data.in.num_lig, data.in.c, typeTemp); // second fib. sphere
        int lig_actual = spherical_wedge(data.in.num_lig, typeTemp, angle, data.in.c);

        Atom patch = Atom(1,1,1,typeLig);
        gen_ligands(data, ligand, patch, typeNano); // find closest spheres on first fib. sphere and change their type

        beads.erase(beads.begin()+data.in.num_of_beads, beads.end()); // erase second fib sphere
    }

protected:
    /**
     * @brief spherical_wedge
     * @param num_beads
     * @param type
     * @param angle - in degrees (0 - 360)
     */
    int spherical_wedge(int num_beads, int type, double angle, double c) {
        int beg = beads.size();
        //
        // generate homogenious distribution
        //
        fibonacci_sphere(beads, num_beads, type);

        //
        // erase beads outside of wedge
        //
        for(int i=beads.size(); i>beg; --i) {
            if(!isWedge(i, angle)) {
                beads.erase(beads.begin()+i);
            }
        }

        //
        // erase beads on poles based on param c
        //
        for(int i=beads.size(); i>beg; --i) {
            if( fabs( beads[i].pos.z ) > (2.0-c)*0.5 ) {
                beads.erase(beads.begin()+i);
            }
        }

        return beads.size() - beg;
    }

    bool isWedge(int i, int angle)
    {
    	Tensor_xyz a_x(1,0,0);
        Tensor_xyz o = beads[i].pos;
        o.z = 0;
        o.normalise();

        if( ( o.cross(a_x) ).size() < sin(angle*degToRad) && o.dot(a_x) > 0.0)
            return true;
        return false;
    }
};

class TennisBall2: public Sphere {
public:
    inline static const string name = "tennis_ball_2";

    TennisBall2() : Sphere() {}

    void generate( Data data )
    {
        Atoms ligand;

        fibonacci_sphere(beads, data.in.num_of_beads, typeNano);
        fibonacci_sphere(ligand, data.in.num_lig, typeTemp); // second fib. sphere
        gen_tennisball(data.in.c, typeTemp); // turn beads on sphere to type to look like tennis ball -> half of first sphere is now type_temp
        Atom patch = Atom(1,1,1,typeLig);
        gen_ligands(data, ligand, patch, typeTemp); // last num_lig beads (second sphere) -> find closest beads on first sphere of type_temp and turn them type_lig

        for(int i=0; i<beads.size(); ++i){
            if(beads[i].type == typeTemp)
                beads[i].type = typeNano;
        }

        beads.erase(beads.begin()+data.in.num_of_beads, beads.end()); // erase second fib sphere
    }

private:
    void gen_tennisball(double S, int type){
        double tolerance = 0.005;
        for(unsigned int i=0; i<beads.size(); ++i) {
            if(tennisball_condition(i,tolerance, S)) {
                beads[i].type = type;
            }
        }
    }

    bool tennisball_condition(unsigned int i, double tolerance, double S) {
        double x0 = cos(S/2.0);
        double y0 = sin(S/2.0);
        double b0 = 1.0/sqrt(2.0);

        //A (x0, y0, 0)
        if( (beads[i].pos.x > x0-tolerance && beads[i].pos.x < x0+tolerance) )
            if( (beads[i].pos.y > y0-tolerance && beads[i].pos.y < y0+tolerance) || (beads[i].pos.y > -y0-tolerance && beads[i].pos.y < -y0+tolerance) )
                if(beads[i].pos.z > 0.0-tolerance && beads[i].pos.z < 0.0+tolerance)
                    return true;

        //B
        if(beads[i].pos.x > 0.0-tolerance && beads[i].pos.x < 0.0+tolerance)
            if( (beads[i].pos.y > b0-tolerance && beads[i].pos.y < b0+tolerance) || (beads[i].pos.y > -b0-tolerance && beads[i].pos.y < -b0+tolerance))
                if( (beads[i].pos.z > b0-tolerance && beads[i].pos.z < b0+tolerance) || (beads[i].pos.z > -b0-tolerance && beads[i].pos.z < -b0+tolerance))
                    return true;

        //C
        if( (beads[i].pos.x > -x0-tolerance && beads[i].pos.x < -x0+tolerance) )
            if(beads[i].pos.y > 0.0-tolerance && beads[i].pos.y < 0.0+tolerance)
                if( (beads[i].pos.z > y0-tolerance && beads[i].pos.z < y0+tolerance) || (beads[i].pos.z > -y0-tolerance && beads[i].pos.z < -y0+tolerance) )
                    return true;

        // function g = connects C to B -> area defined by this curve
        if(beads[i].pos.y < g(beads[i].pos.x, x0, y0, b0) && beads[i].pos.y > -g(beads[i].pos.x, x0, y0, b0)) // Y
            return true;

        // function h = connect B to A -> area not defined by this curve
        if(beads[i].pos.x > 0)
            if(beads[i].pos.y < h(beads[i].pos.x, x0, y0, b0) && beads[i].pos.y > -h(beads[i].pos.x, x0, y0, b0)) // Y
                return true;

        return false;
    }

    double h(double x, double x0, double y0, double b0) {
        return ( x*((y0-b0)/(x0)) ) + b0;
    }

    double g(double x, double x0, double y0, double b0) {
        return sqrt(1.0 - x*x - h(-1.0*x, x0, y0, b0)*h(-1.0*x, x0, y0, b0) );
    }
};

#endif // TENNISBALL_H
