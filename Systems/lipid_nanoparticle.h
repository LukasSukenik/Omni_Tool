#ifndef LIPID_NANOPARTICLE_H
#define LIPID_NANOPARTICLE_H

#include "system_base.h"
#include "atom.h"
#include "xtcanalysis.h"
#include "histogram.h"
#include "cluster_analysis.h"
#include "rdf.h"


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
        ss << help_calc_water_content() << endl;
        ss << help_print_last_frame_as_gro() << endl;
        ss << help_analyze_phase() << endl;
        ss << help_cluster_analysis() << endl;
        return ss.str();
    }

    void execute(Data& data)
    {
        data.in.param.validate_keyword("System_execute", "calc_water_content | print_last_frame_as_gro | analyze_phase | Cluster_Analysis");
        if(data.in.param["System_execute"].compare("calc_water_content") == 0) { calc_water_content(data); }
        if(data.in.param["System_execute"].compare("print_last_frame_as_gro") == 0) { print_last(data); }
        if(data.in.param["System_execute"].compare("analyze_phase") == 0) { analyze_phase(data); }
        if(data.in.param["System_execute"].compare("Cluster_Analysis") == 0) { cluster_analysis(data); }
    }

private:

    ///
    /// Cluster analysis
    ///
    string help_cluster_analysis()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Lipid_Nanoparticle" << endl;
        ss << "System_execute: Cluster_Analysis" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "Atom_type: 2 3 5 6" << endl;
        ss << "Cluster_cutoff: 2.0" << endl;
        ss << "ID: 1" << endl;
        ss << "Optionally use keyword 'Only_last_frame:' , i.e. use only last trajectory frame" << endl;
        ss << "Optionally use keyword 'Trajectory_frame: 10' , i.e. use only 10th frame" << endl;
        ss << "Optionally use keyword 'Trajectory_step: 10' , i.e. only every 10th step of trajectory" << endl;

        return ss.str();
    }

    void validate_cluster_analysis_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "traj_1.xtc");
        data.in.p_vec_int.validate_keyword("Atom_type", "2 3 5 6");
        data.in.p_float.validate_keyword("Cluster_cutoff", "2.0");
    }

    void cluster_analysis(Data& data)
    {
        validate_cluster_analysis_inputs(data);
        cerr << "Lipid_Nanoparticle::cluster_analysis" << endl;

        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& topo = data.coll_beads[sys_id];
        Clusters clusters(topo, data.in.p_vec_int["Atom_type"]); // list of particle indexes
        Trajectory traj(data.in.param["Trajectory_file"]);

        size_t step = (data.in.p_int.contains("Trajectory_step")) ? data.in.p_int["Trajectory_step"] : 1;
        size_t start = (data.in.p_bool.contains("Only_last_frame")) ? traj.frame_count()-1 : 0;

        for(size_t i=start; i<traj.frame_count(); i+=step)
        {
            topo.set_frame(traj[i]);
            clusters.analyze(topo, data.in.sim_box, data.in.p_float["Cluster_cutoff"]);
            cout << i << " " << clusters.size() << " ";
            for(Cluster& cluster : clusters)
            {
                cout << cluster.size()/3 << " "; // dividing by 3 to get the lipid count, assumes we analyze tails only, deserno 4 bead types model
            }
            cout << endl;
            clusters.clear();
        }
    }




    ///
    /// Cluster Radial Distribution Function
    ///
    string help_cluster_rdf()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Lipid_Nanoparticle" << endl;
        ss << "System_execute: Cluster_Analysis" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "Atom_type: 2 3 5 6" << endl;
        ss << "Cluster_cutoff: 2.0" << endl;
        ss << "ID: 1" << endl;
        ss << "Optionally use keyword 'Only_last_frame:'" << endl;
        ss << "Optionally use keyword 'Trajectory_step: 10'" << endl;

        return ss.str();
    }

    void validate_cluster_rdf_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "traj_1.xtc");
        data.in.p_vec_int.validate_keyword("Atom_type", "2 3 5 6");
        data.in.p_float.validate_keyword("Cluster_cutoff", "2.0");
    }

    void cluster_rdf(Data& data)
    {
        validate_cluster_analysis_inputs(data);
        cerr << "Lipid_Nanoparticle::cluster_analysis" << endl;

        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& topo = data.coll_beads[sys_id];
        Atoms clust_topo;
        clust_topo.reserve(topo.size());
        Clusters clusters(topo, data.in.p_vec_int["Atom_type"]); // list of particle indexes
        Trajectory traj(data.in.param["Trajectory_file"]);

        size_t step = (data.in.p_int.contains("Trajectory_step")) ? data.in.p_int["Trajectory_step"] : 1;
        size_t start = (data.in.p_bool.contains("Only_last_frame")) ? traj.frame_count()-1 : 0;
        for(size_t i=start; i<traj.frame_count(); i+=step)
        {
            topo.set_frame(traj[i]);
            clusters.analyze(topo, data.in.sim_box, data.in.p_float["Cluster_cutoff"]);
            cout << i << " " << clusters.size() << " ";
            for(Cluster& cluster : clusters)
            {
                cout << cluster.size()/3 << " "; // dividing by 3 to get the lipid count, assumes we analyze tails only, deserno 4 bead types model
                clust_topo.set_cluster(topo, cluster);
                rdf(clust_topo, 0.0, 50.0, 200, "clust_rdf");
            }
            cout << endl;
            clusters.clear();
        }
    }

    ///
    /// Analyze Phase
    ///
    string help_analyze_phase()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Lipid_Nanoparticle" << endl;
        ss << "System_execute: analyze_phase" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "Histo_2D_dirs_outfile: dirs_distrib_2D" << endl;
        ss << "Histo_1D_dirs_outfile: dirs_distrib_1D" << endl;
        ss << "Histo_spherical_settings: 20 40" << endl;
        ss << "Averaged_frame_count: 10" << endl;
        ss << "ID: 1" << endl;
        ss << "Histo 1D is a cumulative ordered version of histo_2D" << endl;

        return ss.str();
    }

    void validate_analyze_phase_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file",             "data.start");
        data.in.param.validate_keyword("Trajectory_file",       "traj_1.xtc");
        data.in.param.validate_keyword("Histo_2D_dirs_outfile", "dirs_distrib_2D");
        data.in.param.validate_keyword("Histo_1D_dirs_outfile", "dirs_distrib_1D");
        data.in.p_vec_int.validate_keyword("Histo_spherical_settings", "20 40");
        data.in.p_int.validate_keyword("Averaged_frame_count",     "10");
    }

    void analyze_phase(Data& data)
    {
        validate_analyze_phase_inputs(data);
        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& topo = data.coll_beads[sys_id];

        Trajectory traj(data.in.param["Trajectory_file"]);
        string phase = "";

        phase.append( analyze_dirs(data, topo, traj, 10) ); // Dirs analysis identifies planar membrane
        cout << phase << endl;
    }

    ///
    /// Print Last Frame
    ///
    string help_print_last_frame_as_gro()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Lipid_Nanoparticle" << endl;

        return ss.str();
    }

    void print_last(Data& data)
    {
        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& mem = data.coll_beads[sys_id];

        Trajectory traj(data.in.param["Trajectory_file"]);
        mem.set_frame(  traj[traj.frame_count()-1]  );

        data.gro.print_lammps_data(mem, data.in.sim_box.get_box());
    }

    ///
    /// Calc Water Content
    ///
    string help_calc_water_content()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Lipid_Nanoparticle" << endl;

        return ss.str();
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

        Trajectory traj(data.in.param["Trajectory_file"]);

        int sys_id = data.id_map[ data.in.p_int["ID"] ];
        Atoms& last_frame = data.coll_beads[sys_id];
        last_frame.set_frame(traj[ traj.frame_count()-1 ]);
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

    //
    // Other
    //

    void get_dirs(Data& data, Atoms& frame, vector<Tensor_xyz>& dirs)
    {
        Tensor_xyz periodic = Tensor_xyz(data.in.sim_box.xhi - data.in.sim_box.xlo, data.in.sim_box.yhi - data.in.sim_box.ylo, data.in.sim_box.zhi - data.in.sim_box.zlo);
        for(size_t i=0; i<frame.size(); i+=4)
        {
            dirs[i/4] = (frame[i].pos - frame[i+3].pos);

            // fix for periodic box
            if(dirs[i/4].x >  0.5*periodic.x) { dirs[i/4].x -= periodic.x; }
            if(dirs[i/4].x < -0.5*periodic.x) { dirs[i/4].x += periodic.x; }
            if(dirs[i/4].y >  0.5*periodic.y) { dirs[i/4].y -= periodic.y; }
            if(dirs[i/4].y < -0.5*periodic.y) { dirs[i/4].y += periodic.y; }
            if(dirs[i/4].z >  0.5*periodic.z) { dirs[i/4].z -= periodic.z; }
            if(dirs[i/4].z < -0.5*periodic.z) { dirs[i/4].z += periodic.z; }

            dirs[i/4].normalise();
        }
    }

    void get_frame_averaged_dirs(Data& data, vector<Tensor_xyz>& dirs, Atoms& mem, Trajectory& traj, int number_of_averaged_frames)
    {
        vector<Tensor_xyz> dirs_temp(dirs.size(), Tensor_xyz(0.0, 0.0, 0.0));

        for(int i=0; i<number_of_averaged_frames; ++i)
        {
            mem.set_frame(  traj[traj.frame_count()-(1+i)]  );
            get_dirs(data, mem, dirs_temp);
            for(size_t j=0; j<dirs_temp.size(); ++j)
            {
                dirs[j] = dirs[j] + dirs_temp[j];
            }
        }
        for(size_t j=0; j<dirs.size(); ++j)
        {
            dirs[j] = dirs[j]*(1.0/number_of_averaged_frames);
        }
    }

    string analyze_dirs(Data& data, Atoms& topo, Trajectory& traj, int number_of_averaged_frames)
    {
        vector<Tensor_xyz> dirs(topo.size() / 4, Tensor_xyz(0.0, 0.0, 0.0));
        get_frame_averaged_dirs(data, dirs, topo, traj, number_of_averaged_frames);

        Histogram_Spherical h_sp(data.in.p_vec_int["Histo_spherical_settings"][0], data.in.p_vec_int["Histo_spherical_settings"][1]);

        for(Tensor_xyz a : dirs)
        {
            h_sp.add(a, Tensor_xyz(0.0,0.0,1.0), 0.0);
        }
        h_sp.print(data.in.param["Histo_2D_dirs_outfile"]);
        h_sp.print_ordered_cumulative(data.in.param["Histo_1D_dirs_outfile"]);

        if( h_sp.is_planar(dirs.size()) )
        {
            return "planar";
        }
        else
        {
            return "UNK";
        }

        /*Histogram_Spherical h_rot(20,40);
        Tensor_xyz dir_highest = h_sp.get_highest();
        Tensor_xyz axis;
        double angle;
        if(dir_highest.x != 1.0)
        {
            axis = dir_highest.cross(Tensor_xyz(0.0, 1.0, 0.0));
            angle = asin( axis.size() );
            axis.normalise();
            for(Tensor_xyz a : dirs)
            {
                h_rot.add(a, axis, angle*rad_to_deg);
            }
        }
        h_rot.print();
        h_rot.print_ordered();*/
    }
};


#endif // LIPID_NANOPARTICLE_H
