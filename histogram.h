#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <iostream>
#include <fstream>
#include <vector>
#include "tensor_xyz.h"

using namespace std;




class Histogram
{
public:
    vector<int> hist;
    double min, max;

    Histogram(){}

    Histogram(int size, double minn, double maxx)
    {
        this->min = minn;
        this->max = maxx;
        hist = vector<int>(size, 0);
    }

    void add(double value)
    {
        value = value - min;
        int index = (int) ( value * hist.size() / (max-min) );
        hist[index]++;
    }

    size_t size()
    {
        return hist.size();
    }

    void print(string filename)
    {
        fstream fs(filename, fstream::out);
        if(!fs.is_open())
        {
            cout << "Failed to open file " << filename << endl;;
        }
        else
        {
            for(size_t i=0; i<hist.size(); ++i)
            {
                fs << min+ ( ((max-min)*i)/hist.size() )  << " " << hist[i] << endl;
            }
        }
        fs.close();
    }

    double get(std::size_t index)
    {
        return hist[index];
    }
};




class Histogram_2D
{
public:
    vector< vector<int>> hist;
    double min1, min2, max1, max2;

    Histogram_2D(){}

    Histogram_2D(int size1, int size2, double minn1, double maxx1, double minn2, double maxx2)
    {
        this->min1 = minn1;
        this->max1 = maxx1;
        this->min2 = minn2;
        this->max2 = maxx2;
        hist.resize(size1);
        for(size_t i=0; i<hist.size(); ++i)
        {
            hist[i].resize(size2);
            for(size_t j=0; j<hist[i].size(); ++j)
            {
                hist[i][j] = 0;
            }
        }
    }

    void add(double val1, double val2)
    {
        val1 = val1 - min1;
        val2 = val2 - min2;
        int i = (int) ( val1 * hist.size() / (max1-min1) );
        int j = (int) ( val2 * hist[i].size() / (max2-min2) );
        hist[i][j]++;
    }
};




/**
 * @brief The Histogram_Spherical class
 * - histogram of a sphere surface, each slot needs to have equal surface are
 * - separate sphere into sections based on canopy surface and azimuth angle
 * - use canopy surface -> depends on z coords -> split z dir equidistantly
 * - azimuth angle equidistantly
 *
 */
class Histogram_Spherical
{
public:
    const double pi = 3.141592653589793;
    const double deg_to_rad = 0.0174532925199;
    const double rad_to_deg = 57.2957795131;
    int size1;
    int size2;

    Histogram inclination;
    Histogram azimuth;
    Histogram_2D h_2D;

    Histogram_Spherical(int size1, int size2) : size1(size1), size2(size2)
    {
        inclination = Histogram(size1, 0, pi);
        azimuth = Histogram(size2, -pi, pi);
        h_2D = Histogram_2D(size1, size2, -1.0, 1.0, -pi, pi);
    }

    void add(Tensor_xyz dir, Tensor_xyz axis, double angle)
    {
        Tensor_xyz temp = dir;
        if(angle != 0.0)
        {
            temp.rotate(axis, angle*deg_to_rad);
            temp.normalise();
        }
        azimuth.add(calc_phi(temp));
        inclination.add(calc_theta(temp));
        h_2D.add(temp.z,calc_phi(temp));
    }

    Tensor_xyz get_highest()
    {
        vector<int> array = get_ordered_array();
        size_t size_1 = h_2D.hist.size();
        size_t size_2 = h_2D.hist[0].size();
        //double z;
        double rad;
        double phi;
        Tensor_xyz r;

        for(size_t i=0; i<size_1; ++i)
        {
            for(size_t j=0; j<size_2; ++j)
            {
                if(array.back() == h_2D.hist[i][j])
                {
                    r.z = -1.0 + (2.0/size_1)*i;
                    phi = (360.0/size_2)*j;
                }
            }
        }
        rad = sqrt(1.0 - r.z*r.z);
        r.x = rad*cos(phi*deg_to_rad);
        r.y = rad*sin(phi*deg_to_rad);
        r.normalise();
        return r;
    }

    void print(string filename)
    {
        fstream fs(filename, fstream::out);
        if(!fs.is_open())
        {
            cout << "Failed to open file " << filename << endl;;
        }
        else
        {
            if(!h_2D.hist.empty())
            {
                size_t size_1 = h_2D.hist.size();
                size_t size_2 = h_2D.hist[0].size();
                for(size_t i=0; i<size_1; ++i)
                {
                    for(size_t j=0; j<size_2; ++j)
                    {
                        fs << (200.0/size_1)*i -100.0 << " " << (360.0/size_2)*j << " " << h_2D.hist[i][j] << endl;
                    }
                }
            }
        }
        fs.close();
    }

    void print_ordered()
    {
        if(!h_2D.hist.empty())
        {
            fstream o_file("distrib_1D", fstream::out);
            vector<int> array = get_ordered_array();
            for(size_t i=0; i<array.size(); ++i)
            {
                o_file << i << " " << array[i] << endl;
            }
            o_file.close();
        }
    }

    void print_ordered_cumulative(string filename)
    {
        fstream fs(filename, fstream::out);
        if(!fs.is_open())
        {
            cout << "Failed to open file " << filename << endl;;
        }
        else
        {
            if(!h_2D.hist.empty())
            {
                vector<int> array = get_ordered_array();
                int sum=0;
                for(size_t i=0; i<array.size(); ++i)
                {
                    sum += array[i];
                    fs << i << " " << sum << endl;
                }
            }
        }
        fs.close();
    }

    /**
     * @brief get_Normalized_Entropy
     *
     * Reduce a histogram to a single value representing randomness
     * - 1 represents maximum randomness (a continuous Uniform distribution)
     * - 0 represents total certainty (a Dirac delta distribution)
     *
     * H_norm = ( SUM p_i log(p_i) ) / (log K)
     * - K number of histogram bins
     * - p_i = count_i / total_count
     *
     */
    double get_Normalized_Entropy(int total_count)
    {
        if(!h_2D.hist.empty())
        {
            double K = size1*size2;
            double count_i = 0.0;
            double p_i = 0.0;
            double H_norm = 0.0;

            size_t size_1 = h_2D.hist.size();
            size_t size_2 = h_2D.hist[0].size();
            for(size_t i=0; i<size_1; ++i)
            {
                for(size_t j=0; j<size_2; ++j)
                {
                    if(h_2D.hist[i][j] > 0)
                    {
                        count_i = h_2D.hist[i][j];
                        p_i = count_i / total_count;
                        H_norm -= p_i * log(p_i) / log(K);
                    }
                }
            }
            return H_norm;
        }
        return -1.0;
    }

    void test_uniformity()
    {
        Tensor_xyz temp;
        for(int i=0; i<100*1000*1000; ++i)
        {
            temp.randomUnitSphere();
            this->add(temp, Tensor_xyz(0.0, 0.0, 1.0), 0.0);
        }
    }

private:
    double calc_phi(Tensor_xyz& dir)
    {
        double phi = 0.0; // -pi to pi
        if(dir.x >  0.0) phi = atan(dir.y/dir.x);
        if(dir.x <  0.0 && dir.y >= 0.0) phi = atan(dir.y/dir.x) + pi;
        if(dir.x <  0.0 && dir.y <  0.0) phi = atan(dir.y/dir.x) - pi;
        if(dir.x == 0.0 && dir.y >  0.0) phi =  pi * 0.5;
        if(dir.x == 0.0 && dir.y <  0.0) phi = -pi * 0.5;
        if(dir.x == 0.0 && dir.y == 0.0 && (dir.z == 1.0 || dir.z == -1.0) ) phi = 0.0;
        if(dir.x == 0.0 && dir.y == 0.0 && !(dir.z == 1.0 || dir.z == -1.0) )
        {
            cerr << "Histogram_Spherical::phi angle undefined -> " << dir << endl;
            exit(-1);
        }
        return phi;
    }

    double calc_theta(Tensor_xyz& dir)
    {
        double theta = 0.0; // 0.0 to pi
        if(dir.z > 0.0) theta = atan( sqrt(dir.x*dir.x + dir.y*dir.y) / dir.z );
        if(dir.z < 0.0) theta = pi + atan( sqrt(dir.x*dir.x + dir.y*dir.y) / dir.z );
        if(dir.z == 0.0 && sqrt(dir.x*dir.x + dir.y*dir.y) != 0.0 ) theta = 0.5*pi;
        if(dir.x == 0.0 && dir.y == 0.0 && dir.z == 0.0 )
        {
            cerr << "Histogram_Spherical::theta angle undefined -> " << dir << endl;
            exit(-1);
        }
        return theta;
    }

    vector<int> get_ordered_array()
    {
        size_t size_1 = h_2D.hist.size();
        size_t size_2 = h_2D.hist[0].size();
        vector<int> array(size_1 * size_2, 0);
        for(size_t i=0; i<size_1; ++i)
        {
            for(size_t j=0; j<size_2; ++j)
            {
                array[i*size_2 + j] = h_2D.hist[i][j];
            }
        }
        sort(array.begin(), array.end()); // default sort operator <

        return array;
    }
};

#endif // HISTOGRAM_H
