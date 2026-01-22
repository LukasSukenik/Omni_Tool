#ifndef OBLATESPHEROID_H
#define OBLATESPHEROID_H

#include <complex>

#include "sphere.h"

using namespace std;

/**
 * @brief The OblateSpheroid class - generate oblate spheroid with surface area of 4*PI (equivalent to sphere of unit radius)
 */
class OblateSpheroid: public Sphere {
public:
    inline static const string keyword = "oblate_spheroid";
    const string name = "oblate_spheroid";

    OblateSpheroid() : Sphere("oblate_spheroid") {}

    /**
     * @brief generate (x^2+y^2)/(a^2) + (z^2/c^2) = 1, c > 1 prolate, a>1 oblate
     * @param size
     */
    void generate( Data& data )
    {
        vector<Atom> ligand;

        typeNano = data.in.atom_type;
        typeLig = data.in.atom_type + 1;

        int orientations = orientY;
        int num_beads = data.in.beads_per_area * surface(data.in.scale, data.in.scale*data.in.c);
        int num_ligs = data.in.ligs_per_area * surface(data.in.scale, data.in.scale*data.in.c);

        fibonacci_spheroid(beads, num_beads, data.in.c, typeNano, orientations);
        fibonacci_spheroid(ligand, num_ligs, data.in.c, typeTemp, orientations);
        gen_ligands( data, ligand, data.in.patch_1, typeNano, typeTemp);
        gen_ligands( data, ligand, data.in.patch_2, typeNano, typeTemp);

        gen_TMD(data.in.c, data.in.scale, data.in.tmd.size, data.in.tmd.proximal_n, data.in.tmd.distal_n, orientations);

        // add to existing data
        int i=0;
        for(auto& item : beads)
        {
            item.mol_tag = data.in.mol_tag;
            item.N = i+1+data.in.offset+data.get_bead_count();
            ++i;
        }
        i=0;
        for(auto& b : bonds)
        {
            b.N += data.get_bond_count();
            b.at1 += data.get_bead_count();
            b.at2 += data.get_bead_count();
        }
        for(auto& a : angles)
        {
            a.at1 += data.get_bead_count();
            a.at2 += data.get_bead_count();
            a.at3 += data.get_bead_count();
        }

        cerr << "Oblate beads size " << beads.size() << endl;
    }



    double surface(double equator_r, double polar_r)
    {
        const double PI = 3.141592653589793;
        double a = equator_r;
        double c = polar_r;
        if(polar_r > equator_r)
        {
            double frac = (a*a) / (c*c);
            double e = sqrt(1 - frac );
            return 2*PI*a*a * (1+ c/(a*e)*asin(e));
        }
        return 0.0;
    }

    /**
     * @brief gen_TMD - generate transmembrane domain
     */
    void gen_TMD(double c, double scale, int size, int n_proxima, int n_distal, int orientation)
    {
        double increment = 0.9 * 1.0/scale;
        if(orientY == orientation)
        {
            int type;
            beads.push_back( Atom(0, -c+increment, 0, 8) );
            beads.push_back( Atom(0, -c, 0, 8) );
            bonds.push_back( Bond(1, 6, beads.size()-1, beads.size() ) ); // beads indexed from 0, but the at1, at2, values indexed from 1 -> so no -1

            for(int i=1; i<=size; ++i)
            {
                type = 12;
                if(i <= n_proxima || i > size-n_distal) { type = 11; }
                beads.push_back( Atom(0, -c -i*increment, 0, type) );
                bonds.push_back(Bond(i+1, 6, beads.size()-1, beads.size() ));
                angles.push_back( Angle(i, 1, beads.size()-2, beads.size()-1, beads.size() ) );
            }
        }
    }

protected:
    void fibonacci_spheroid(vector<Atom>& container, int samples, double c, const int type, int orientation)
    {
        const double PI = 3.141592653589793;
        double tolerance = 0.0001;

        double offset = c*2.0/samples;
        double increment = PI * ( 3.0 - sqrt(5.0));

        double x,y,z,r,phi;
        double temp;

        for(int i=0; i<samples; ++i) {

            z = ((i * offset) - c) + (offset * 0.5); // z linearly distributed from -c to c

            //
            // c ~ 0.0 -> disk, badly described by spheroid function
            // c ~ 1.0 -> sphere, cant describe by disk function
            //
            // Here we scale between the two function based on c parameter
            //
            temp = (1.0-c) * ( c*sqrt(fabs(z))/sqrt(c) ); // disk equation
            if(z < 0.0)
                temp *= -1.0;
            temp += c * transform2(z,c);

            //
            // we further transform the function to purely spheroid function for nanoparticle sides, which are again badly described by disk function
            //
            if( c > 1.0 || fabs( transform2(z,c) ) < c*sqrt(fabs(z))/sqrt(c) )
                temp = transform2(z,c);

            //cout << z << " " << temp << endl;
            z = temp;

            // prescribe spheroid
            r = sqrt( 1 - ((z*z)/(c*c)) ); // r = sqrt(1 - z*z); -> c=1.0 -> sqrt(1-z*z) - z*z
            phi = i * increment;

            x = cos(phi) * r;
            y = sin(phi) * r;

            if( ( ((x*x)+(y*y)) + ((z*z)/(c*c)) ) > (1 - tolerance) && ( ((x*x)+(y*y)) + ((z*z)/(c*c)) ) < (1 + tolerance)) {
                if(orientX == orientation)
                    container.push_back(Atom(z,x,y,type));
                if(orientY == orientation)
                    container.push_back(Atom(x,z,y,type));
                if(orientZ == orientation)
                    container.push_back(Atom(x,y,z,type));
            }
        }
        //exit(0);
    }

private:
    //
    // STEPS TO CREATE WHATEVER DISTRIBUTION:
    //
    // DECIDE on probability density, integrate probability density -> get cumulative distribution function, inverse the function
    //

    //
    // r(z) = sqrt(1-(z*z)/(c*c)) -> radius as a function of z
    //

    //
    // Oblate spheroid: S = 2 * pi * integral{ r(z) * sqrt[ 1 + derivate( r(z) )^2 ] } dz
    //      derivate r * sqrt(1+ (derivate r)^2 ) = sqrt(  1 + (1-c*c)*z*z/(c*c*c*c)  )
    //

    //
    // Cumulative distribution function: 2 * pi * { integrate sqrt(  1 + (1-c^2)*z^2 / (c^4) ) dz }
    //
    // Cumulative distribution function: 2 * pi * { integrate     (  1 + (1-c*c)*z^2 / (c^4) )^(0.5) dz }
    //
    // pi * (z * sqrt(  1 + z*z/(c*c*c*c) - z*z/(c*c)  )   +   (c*c * sinh_inv(  ( sqrt(1-c*c) * z) / (c*c) )  ) / sqrt(1 - c*c))
    //
    // pi * (z * sqrt(  1 + z*z/(c*c*c*c) - z*z/(c*c)  )   +   (c*c * sinh^-1 (  ( sqrt(1-c*c) * z) / (c*c) )  ) / sqrt(1 - c*c))
    //
    // + inverse entire function -> mirror function on function y=x
    //

    /**
     * @brief cumulDistFce Cumulative distribution function of proplate/oblate spheroid for a=1 in z dir
     *
     * http://mathworld.wolfram.com/ProlateSpheroid.html
     * http://mathworld.wolfram.com/OblateSpheroid.html
     *
     * @param z
     * @param c
     * @return
     */
    std::complex<double> cumulDistFce(std::complex<double> z, std::complex<double> c)
    {
        return (  z * sqrt(1.0 + z*z/(c*c*c*c) - z*z/(c*c))  +  ( c*c * asinh(( sqrt(1.0-c*c)*z ) / (c*c)) ) / sqrt(1.0-c*c)  ) / norm(c);
    }

    /**
     * @brief norm
     * @param c - parameter for norming function from proplate/oblate spheroid surface term
     * @return
     */
    std::complex<double> norm(std::complex<double> c)
    {
        return (1.0 + (c*c * asinh( sqrt(1.0-c*c)/c ) ) / sqrt( 1.0 - c*c ) ) / c;
    }

    /*double transform(double z, double c) {
        return (z > 0.0) ? ( cumulDistFce(z-c,c) + c ) : ( cumulDistFce(z+c,c) - c );
    }*/

    double transform2(double z, double c)
    {
        std::complex<double> zz(z,0.0);
        std::complex<double> cc(c,0.0);

        return real(inverse(cumulDistFce(zz,cc), zz));
    }

    std::complex<double> inverse(std::complex<double> value, std::complex<double> xy)
    {
        return xy + xy - value;
    }
};

#endif // OBLATESPHEROID_H
