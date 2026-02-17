#ifndef XTCANALYSIS_H
#define XTCANALYSIS_H

#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <vector>

#include "tensor_xyz.h"

#include "xdrfile-1.1.4/include/xdrfile.h"
#include "xdrfile-1.1.4/include/xdrfile_xtc.h"

using namespace std;

class Trajectory : public vector< vector<Tensor_xyz> >
{
private:
    char fileName[64];
    int natoms;

    int step;
    int release_step = numeric_limits<int>::max();
    float time;
    matrix temp_box;
    float prec;

    int first_step;
    int status=exdrOK;

public:   
    vector<Tensor_xyz> box_traj;

    Trajectory() {}

    void load(string inName, int start=-1, int stop=-1)
    {
        Tensor_xyz pos;

        string_to_char(inName); // inName -> fileName

        status = read_xtc_natoms(fileName, &natoms);
        if (status == exdrOK)
        {
            XDRFILE* xfp = xdrfile_open(fileName, "r");
            if (xfp != NULL)
            {
                rvec k[natoms];

                status = read_xtc(xfp, natoms, &first_step, &time, temp_box, k, &prec);

                box_traj.push_back(Tensor_xyz(temp_box[0][0], temp_box[1][1], temp_box[2][2]));
                this->push_back( vector<Tensor_xyz>() );
                for(int i=0; i<natoms; ++i)
                {
                    this->back().resize(natoms);
                    this->back()[i].x = k[i][0];
                    this->back()[i].y = k[i][1];
                    this->back()[i].z = k[i][2];
                }

                while(status == exdrOK && ( stop == -1 || step < first_step+stop) ) // Load frame
                {
                    status = read_xtc(xfp, natoms, &step, &time, temp_box, k, &prec);

                    if(step > start)// add frame
                    {
                        box_traj.push_back(Tensor_xyz(temp_box[0][0], temp_box[1][1], temp_box[2][2]));
                        this->push_back( vector<Tensor_xyz>() );
                        for(int i=0; i<natoms; ++i)
                        {
                            this->back().resize(natoms);
                            this->back()[i].x = k[i][0];
                            this->back()[i].y = k[i][1];
                            this->back()[i].z = k[i][2];
                        }
                    }
                }
                xdrfile_close(xfp);
            }
            else
            {
                cout << "XTC_Loader::load - File not opened:" << fileName << endl;
            }
        }
        else
        {
            cout << "XTC_Loader::load read_xtc_natoms failure; return code " << status << endl;
            cout << "Filename: " << inName << endl;
        }
    }

    inline int frame_count()
    {
        return size();
    }

private:

    void string_to_char(string inName)
    {
        int i;
        for (i = 0; inName[i] != '\0'; ++i) {
          fileName[i] = inName[i];
        }
        fileName[i] = '\0';
    }
};

#endif // XTCANALYSIS_H
