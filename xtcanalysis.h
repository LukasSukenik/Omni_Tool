#ifndef XTCANALYSIS_H
#define XTCANALYSIS_H

#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <vector>

#include "tensor_xyz.h"
#include "data.h"

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
    int status=exdrOK;

    bool only_last=false;

    int release_step = numeric_limits<int>::max();

public:   
    vector< vector<Tensor_xyz> > conf_traj;
    vector<float> time_traj;
    vector<int > step_traj;
    vector<Tensor_xyz> box_traj;

    Trajectory(){}
    Trajectory(Data& data)
    {
        int start=0;
        int modulo=1;
        int stop=1000*1000;

        if( data.in.p_vec_int.contains("Trajectory_settings") )
        {
            start=data.in.p_vec_int["Trajectory_settings"][0];
            modulo=data.in.p_vec_int["Trajectory_settings"][1];
            stop=data.in.p_vec_int["Trajectory_settings"][2];
        }

        if(data.in.p_bool.contains("Only_last_frame") )
        {
            only_last=true;
        }

        load(data.in.param["Trajectory_file"], start, stop, modulo);
    }
    Trajectory(string inName, int start=0, int stop=1000000, int modulo=1)
    {
        load(inName, start, stop, modulo);
    }

    void load(string inName, int start=0, int stop=1000000, int modulo=1)
    {
        string_to_char(inName); // inName -> fileName
        status = read_xtc_natoms(fileName, &natoms);
        if (status == exdrOK)
        {
            XDRFILE* xfp = xdrfile_open(fileName, "r");
            if (xfp != NULL)
            {
                rvec conf[natoms];
                int frame_count = 0; // for start, stop, modulo
                while(status == exdrOK && (frame_count<1000*1000) ) // Load frame
                {
                    status = read_xtc(xfp, natoms, &step, &time, temp_box, conf, &prec);

                    if(only_last)
                    {
                        save_last(conf);
                    }
                    else
                    {
                        save(conf); // step, time, temp_box are global
                    }
                    ++frame_count;
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

    inline int get_step(int i)
    {
        return step_traj[i];
    }

    inline size_t frame_count()
    {
        return conf_traj.size();
    }

private:


    void save(rvec *conf)
    {
        time_traj.push_back(time);
        step_traj.push_back(step);
        box_traj.push_back(Tensor_xyz(temp_box[0][0], temp_box[1][1], temp_box[2][2])); // box
        conf_traj.push_back( vector<Tensor_xyz>(natoms) );
        for(int i=0; i<natoms; ++i) // conf
        {
            conf_traj.back()[i].x = conf[i][0];
            conf_traj.back()[i].y = conf[i][1];
            conf_traj.back()[i].z = conf[i][2];
        }
    }

    void save_last(rvec *conf)
    {
        if(conf_traj.empty())
        {
            save(conf);
        }
        else
        {
            overwrite(conf);
        }
    }

    void overwrite(rvec *conf)
    {
        time_traj.back() = time;
        step_traj.back() = step;
        box_traj.back() = Tensor_xyz(temp_box[0][0], temp_box[1][1], temp_box[2][2]); // box
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
