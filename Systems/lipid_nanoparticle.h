#ifndef LIPID_NANOPARTICLE_H
#define LIPID_NANOPARTICLE_H

#include "system_base.h"
#include "atom.h"
#include "xtcanalysis.h"

class Histogram
{
public:
    vector<int> hist;
    double min, max;

    Histogram(int size, double minn, double maxx)
    {
        this->min = minn;
        this->max = maxx;
        hist = vector<int>(size, 0);
    }

    void add(double value)
    {
        value = value - min;
        int index = (int) ( value * hist.size() / (max-min) + 1e-9 ); // max would be out-of-bound
        hist[index]++;
    }

    int size()
    {
        return hist.size();
    }
};

class Histogram_Radial : public Histogram
{
public:
    Histogram_Radial(int size) : Histogram(size, -3.141592653589793, 3.141592653589793) {}

    void add(Tensor_xyz dir, Tensor_xyz axis)
    {
        double theta = 0.0; // -pi to pi
        if(dir.x >  0.0) theta = atan(dir.y/dir.x);
        if(dir.x <  0.0 && dir.y >= 0.0) theta = atan(dir.y/dir.x) + 3.141592653589793;
        if(dir.x <  0.0 && dir.y <  0.0) theta = atan(dir.y/dir.x) - 3.141592653589793;
        if(dir.x == 0.0 && dir.y >  0.0) theta =  3.141592653589793 * 0.5;
        if(dir.x == 0.0 && dir.y <  0.0) theta = -3.141592653589793 * 0.5;
        if(dir.x == 0.0 && dir.y == 0.0)
        {
            cerr << "Histogram_Radial::Impossible" << endl;
            exit(-1);
        }
        Histogram::add(theta);
    }
};

class Lipid_Nanoparticle : public System_Base
{
public:
    inline static const string keyword = "Lipid_Nanoparticle";
    const string name = "Lipid_Nanoparticle";

    Force_Field ff;

    Lipid_Nanoparticle() : System_Base("Lipid_Nanoparticle") {}

    string help()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Lipid_Nanoparticle" << endl;

        return ss.str();
    }

    void execute(Data& data)
    {
        if(data.in.system_function.compare("calc_water_content") == 0) { calc_water_content(data); }
        if(data.in.system_function.compare("print_last_frame_as_gro") == 0) { print_last(data); }
        if(data.in.system_function.compare("analyze_phase") == 0) { analyze_phase(data); }
    }

    void analyze_phase(Data& data)
    {
        int sys_id = data.id_map[ data.in.param_int["ID"] ];
        Atoms& mem = data.coll_beads[sys_id];

        Trajectory traj;
        traj.load(data.in.param["Trajectory_file"]);
        mem.set_frame(  traj[traj.frame_count()-1]  );

        vector<Tensor_xyz> dirs;
        dirs.resize(mem.size() / 4);

        for(int i=0; i<mem.size(); i+=4)
        {
            dirs[i/4] = (mem[i].pos - mem[i+3].pos);
            dirs[i/4].normalise();
        }

        Histogram_Radial h_x(360);
        Histogram_Radial h_y(360);
        Histogram_Radial h_z(360);

        for(Tensor_xyz a : dirs)
        {
            h_x.add(Tensor_xyz(a.z, a.y, a.x), Tensor_xyz(1.0, 0.0, 0.0));
            h_y.add(Tensor_xyz(a.x, a.z, a.y), Tensor_xyz(0.0, 1.0, 0.0));
            h_z.add(Tensor_xyz(a.x, a.y, a.z), Tensor_xyz(0.0, 0.0, 1.0));
        }

        for(int i=0; i<h_x.size(); ++i)
        {
            cout << i << " " << h_x.hist[i] << " " << h_y.hist[i] << " " << h_z.hist[i] << endl;
        }


        exit(1);
    }

    void print_last(Data& data)
    {
        int sys_id = data.id_map[ data.in.param_int["ID"] ];
        Atoms& mem = data.coll_beads[sys_id];

        Trajectory traj;
        traj.load(data.in.param["Trajectory_file"]);
        mem.set_frame(  traj[traj.frame_count()-1]  );

        data.gro.print_lammps_data(mem, data.in.sim_box.get_box());
    }

    void calc_water_content(Data& data)
    {
        cerr << "calc_water_content::begin" << endl;

        // feeler atom
        Atom feeler;
        feeler.pos = Tensor_xyz(0.0, 0.0, 0.0);
        feeler.type = 2; // lipid tail

        // set ff
        ff.lj[1] = LJ(-0.8);
        ff.lj[2] = LJ(2.8);
        ff.lj[3] = LJ(2.8);

        // box boundaries
        cerr << data.in.sim_box.xlo << " " << data.in.sim_box.xhi << endl;

        // trajectory load
        Trajectory traj;
        traj.load(data.in.param["Trajectory_file"]);

        int sys_id = data.id_map[ data.in.param_int["ID"] ];
        Atoms& last_frame = data.coll_beads[sys_id];
        last_frame.set_frame(traj[ traj.size()-1 ]);
        cerr << "Trajectory loaded" << endl;

        double step_size = 1.0;
        int water_count = 0;
        int membrane_count = 0;

        double range = data.in.sim_box.xhi - data.in.sim_box.xlo;
        for(double x=data.in.sim_box.xlo; x<data.in.sim_box.xhi; x+= step_size) {
            cerr << "\r" << 100.0 * (x - data.in.sim_box.xlo) / range << " %         " << flush;
            for (double y=data.in.sim_box.ylo; y<data.in.sim_box.yhi; y+=step_size) {
                for (double z=data.in.sim_box.zlo; z<data.in.sim_box.zhi; z+= step_size) {
                    feeler.pos = Tensor_xyz(x, y, z);
                    if( is_water(last_frame, feeler) )
                    {
                        water_count++;
                    }
                    else
                    {
                        membrane_count++;
                    }
                }
            }
        }

        if (membrane_count != 0) {
            double ratio = static_cast<double>(water_count) / membrane_count;
            cout << "Water to Membrane ratio:" << ratio << endl;
            cout << "Membrane %:" << 100.0 * membrane_count / (membrane_count + water_count) << endl;
            cout << "step size:" << step_size << endl;
        } else {
            cerr << "Membrane count is zero, cannot compute ratio." << endl;
        }

        cerr << "calc_water_content::end" << endl;
    }

    bool is_water(Atoms& frame, Atom& feeler)
    {
        return !frame.is_overlap(feeler, 1.0);
    }
};


#endif // LIPID_NANOPARTICLE_H
