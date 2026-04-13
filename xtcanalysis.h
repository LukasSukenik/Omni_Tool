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

class Trajectory
{
private:
    char fileName[64];
    int natoms;

    int step;
    float time;
    matrix temp_box;
    float prec;
    int first_step=-1;
    int status=exdrOK;

    int release_step = numeric_limits<int>::max();

public:   
    vector< vector<Tensor_xyz> > conf_traj;
    vector<float> time_traj;
    vector<int > step_traj;
    vector<Tensor_xyz> box_traj;

    Trajectory(){}
    Trajectory(string inName, int start=-1, int stop=-1, int modulo=1)
    {
        load(inName, start, stop, modulo);
    }

    void load(string inName, int start=-1, int stop=-1, int modulo=1)
    {
        string_to_char(inName); // inName -> fileName
        status = read_xtc_natoms(fileName, &natoms);
        if (status == exdrOK)
        {
            XDRFILE* xfp = xdrfile_open(fileName, "r");
            if (xfp != NULL)
            {
                rvec conf[natoms];
                while(status == exdrOK && ( stop == -1 || step < first_step+stop) ) // Load frame
                {
                    status = read_xtc(xfp, natoms, &step, &time, temp_box, conf, &prec);

                    if(step > start)// add frame
                    {
                        save(conf);
                    }
                }
                xdrfile_close(xfp);
            }
            else
            {
                cout << "Trajectory::load, File not opened: " << fileName << endl;
            }
        }
        else
        {
            if(status == exdrFILENOTFOUND)
            {
                cout << "Trajectory::load, file " << inName << " not found" << endl;
            }
            else
            {
                cout << "Trajectory::load, file " << inName << " MASSIVE, yet mysterious, ERROR" << endl;
            }
            exit(-1);
        }
    }

    vector<Tensor_xyz>& operator[](std::size_t index)
    {
            return conf_traj[index];
    }

    const vector<Tensor_xyz>& operator[](std::size_t index) const
    {
            return conf_traj[index];
    }

    inline size_t frame_count()
    {
        return conf_traj.size();
    }

private:


    void save(rvec *conf)
    {
        box_traj.push_back(Tensor_xyz(temp_box[0][0], temp_box[1][1], temp_box[2][2])); // box
        conf_traj.push_back( vector<Tensor_xyz>(natoms) );
        for(int i=0; i<natoms; ++i) // conf
        {
            conf_traj.back()[i].x = conf[i][0];
            conf_traj.back()[i].y = conf[i][1];
            conf_traj.back()[i].z = conf[i][2];
        }
    }

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
