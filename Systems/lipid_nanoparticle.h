#ifndef LIPID_NANOPARTICLE_H
#define LIPID_NANOPARTICLE_H

#include <algorithm>

#include "system_base.h"
#include "atom.h"
#include "xtcanalysis.h"
#include "histogram.h"
#include "cluster_analysis.h"
#include "rdf.h"

using namespace std;


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
        ss << help_cluster_rdf() << endl;
        ss << help_detect_bleb() << endl;
        ss << help_percolation() << endl;

        ss << help_optional() << endl;

        return ss.str();
    }

    void execute(Data& data)
    {
        data.in.param.validate_keyword("System_execute", "calc_water_content | print_last_frame_as_gro | analyze_phase | Cluster_Analysis");
        if(data.in.param["System_execute"].compare("calc_water_content") == 0) { calc_water_content(data); }
        if(data.in.param["System_execute"].compare("print_last_frame_as_gro") == 0) { print_last(data); }

        if(data.in.param["System_execute"].compare("analyze_phase") == 0) { analyze_phase(data); }
        if(data.in.param["System_execute"].compare("Cluster_Analysis") == 0) { cluster_analysis(data); }
        if(data.in.param["System_execute"].compare("Cluster_RDF") == 0) { cluster_rdf(data); }
        if(data.in.param["System_execute"].compare("Cluster_Surf") == 0) { cluster_surf(data); }
        if(data.in.param["System_execute"].compare("Detect_Bleb") == 0) { detect_bleb(data); }
        if(data.in.param["System_execute"].compare("Percolation_Dim") == 0) { percolation(data); }
    }

private:

    string help_optional()
    {
        stringstream ss;

        ss << "Optionally use keyword 'Only_last_frame:' , i.e. use only last trajectory frame" << endl;
        ss << "Optionally use keyword 'Trajectory_settings: start step stop' , analyze custom frames, like seq" << endl;

        return ss.str();
    }

    ///
    /// Detect Bleb
    ///
    string help_detect_bleb()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Lipid_Nanoparticle" << endl;
        ss << "System_execute: Detect_Bleb" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_detect_bleb_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "traj_1.xtc");
    }

    void detect_bleb(Data& data)
    {
        cerr << "Lipid_Nanoparticle::detect_bleb" << endl;
        validate_detect_bleb_inputs(data);

        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];
        Trajectory traj(data);

        for(size_t ii=0; ii<traj.frame_count(); ++ii)
        {
            topo.set_frame(traj[ii]);
            cout << traj.step_traj[ii] << " " << get_bilayer_lipid_head_set(data, topo).size() << endl;
        }
    }

    ///
    /// Percolation dimensionality
    ///
    string help_percolation()
    {
        stringstream ss;

        ss << "*********************************************************" << endl;
        ss << "System_type: Lipid_Nanoparticle" << endl;
        ss << "System_execute: Percolation_Dim" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "Only_last_frame:" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_percolation_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file",             "data.start");
        data.in.param.validate_keyword("Trajectory_file",       "traj_1.xtc");
    }

    void percolation(Data& data)
    {
        cerr << "Lipid_Nanoparticle::percolation" << endl;
        validate_percolation_inputs(data);

        Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        Trajectory traj(data);
        Lattice_3D grid_3D(data, 1.2);

        for(size_t ii=0; ii<traj.frame_count(); ++ii) // loop over trajectory
        {
            topo.set_frame(traj[ii]);
            grid_3D.add_particles(topo);
            cout << "percolation_dim " << grid_3D.is_percolation( grid_3D.set_start_plane( Tensor_xyz_int(1,0,0) ), Tensor_xyz_int(1,0,0) )+
                                              grid_3D.is_percolation( grid_3D.set_start_plane( Tensor_xyz_int(0,1,0) ), Tensor_xyz_int(0,1,0) )+
                                              grid_3D.is_percolation( grid_3D.set_start_plane( Tensor_xyz_int(0,0,1) ), Tensor_xyz_int(0,0,1) ) << endl;
        }
    }

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
        cerr << "Lipid_Nanoparticle::cluster_analysis" << endl;
        validate_cluster_analysis_inputs(data);

        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];
        Trajectory traj(data);


        for(size_t i=0; i<traj.frame_count(); ++i)
        {
            topo.set_frame(traj[i]);
            Clusters clusters(topo, data.in.p_vec_int["Atom_type"], data.in.sim_box, data.in.p_float["Cluster_cutoff"]); // list of particle indexes, constructor is suepr cheap
            cout << i << " " << traj.get_step(i) << " " << clusters << endl;
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
        ss << "System_execute: Cluster_RDF" << endl;
        ss << "Input_type: lammps_full" << endl;
        ss << "Load_file: data.start" << endl;
        ss << "Trajectory_file: traj_1.xtc" << endl;
        ss << "RDF_outfile: clust_rdf" << endl;
        ss << "Atom_type: 2 3 5 6" << endl;
        ss << "Cluster_cutoff: 2.0" << endl;
        ss << "Only_last_frame: true" << endl;
        ss << "ID: 1" << endl;

        return ss.str();
    }

    void validate_cluster_rdf_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "traj_1.xtc");
        data.in.param.validate_keyword("RDF_outfile", "clust_rdf");
        data.in.p_vec_int.validate_keyword("Atom_type", "2 3 5 6");
        data.in.p_float.validate_keyword("Cluster_cutoff", "2.0");
    }

    void cluster_rdf(Data& data)
    {
        cerr << "Lipid_Nanoparticle::cluster_rdf" << endl;
        validate_cluster_rdf_inputs(data);

        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];
        Trajectory traj(data);
        Clusters clusters(topo, data.in.p_vec_int["Atom_type"]); // list of particle indexes
        RDF rdf(0.0, 60.0, 30);

        Atoms clust_topo;
        Atoms mol_1;
        Atoms mol_2;


        for(size_t i=0; i<traj.frame_count(); ++i)
        {
            topo.set_frame(traj[i]);
            clusters.analyze(topo, data.in.sim_box, data.in.p_float["Cluster_cutoff"]);
            cout << i << " " << traj.get_step(i) << " " << clusters.size() << " ";

            for(Cluster& cluster : clusters)
            {
                cout << cluster.size()/3 << " "; // dividing by 3 to get the lipid count, assumes we analyze tails only, deserno 4 bead types model
                restore_full_lipids(topo, cluster, clust_topo);
                // TODO move CM to 0.0,0.0,0.0, PBC: make cluster whole
                mol_1 = clust_topo.get_molecule(1);
                mol_2 = clust_topo.get_molecule(2);

                rdf.calc(clust_topo);
                rdf.calc(mol_1);
                rdf.calc(mol_2);
                rdf.print(data.in.param["RDF_outfile"]);
            }

            cout << endl;
            clusters.clear();
        }
    }

    ///
    /// Cluster Surface analysis
    /// -> % of Ionizable vs Helper lipids
    /// -> Mixed or Raft
    ///
    string help_cluster_surf()
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

        return ss.str();
    }

    void validate_cluster_surf_inputs( Data& data )
    {
        data.in.param.validate_keyword("Load_file", "data.start");
        data.in.param.validate_keyword("Trajectory_file", "traj_1.xtc");
        data.in.p_vec_int.validate_keyword("Atom_type", "2 3 5 6");
        data.in.p_float.validate_keyword("Cluster_cutoff", "2.0");
    }

    void cluster_surf(Data& data)
    {
        cerr << "Lipid_Nanoparticle::cluster_rdf" << endl;
        validate_cluster_surf_inputs(data);

        Atoms& topo = data.coll_beads[  data.id_map[ data.in.p_int["ID"] ]  ];
        Clusters clusters(topo, data.in.p_vec_int["Atom_type"]); // list of particle indexes
        Trajectory traj(data);
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
        cerr << "Lipid_Nanoparticle::analyze_phase" << endl;
        validate_analyze_phase_inputs(data);

        Atoms& topo = data.coll_beads[   data.id_map[ data.in.p_int["ID"] ]   ];
        Trajectory traj(data);

        //
        // lipid direction randomness - 1 max randomness, 0 total certainty (a Dirac delta distribution)
        //
        /*Histogram_Spherical h_sp = analyze_dirs(data, topo, traj, data.in.p_int["Averaged_frame_count"]);
        h_sp.print(data.in.param["Histo_2D_dirs_outfile"]);
        h_sp.print_ordered_cumulative(data.in.param["Histo_1D_dirs_outfile"]);
        cout << "Enthropy " << h_sp.get_Normalized_Entropy(topo.size()/4);*/

        // Periodic dimensionality: 0-3

        // Amount of Water pockets

        cout << endl;
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
    void restore_full_lipids(Atoms& topo, vector<int> cluster, Atoms& clust_topo)
    {
        clust_topo.resize(cluster.size() * 4 / 3, Atom( Tensor_xyz(0.0,0.0,0.0), 0) );

        sort(cluster.begin(), cluster.end());

        size_t h_c = 0; // lipid head count
        for(size_t i=0; i<cluster.size(); ++i) // loop over tails in cluster
        {
            add_lipid_head(topo, clust_topo, i, h_c, cluster[i]-1);
            add_lipid_middle_bead(topo, clust_topo, i, h_c, cluster[i]);
            add_lipid_last_bead(topo, clust_topo, i, h_c, cluster[i]);
        }
    }

    void add_lipid_head(Atoms& topo, Atoms& clust_topo, size_t& i, size_t& h_c, int p_topo_ID)
    {
        if( (topo[p_topo_ID].type == 1 || topo[p_topo_ID].type == 4) && p_topo_ID >= 0 && i+h_c < clust_topo.size())
        {
            //cerr << "Head [" << i+h_c << "] == [" << p_topo_ID << "] type 1|4==" << topo[ p_topo_ID ].type << endl;
            clust_topo[i+h_c] = topo[ p_topo_ID ]; // write head
            h_c++;
        }
    }

    void add_lipid_middle_bead(Atoms& topo, Atoms& clust_topo, size_t& i, size_t& h_c, int p_topo_ID)
    {
        if( (topo[p_topo_ID].type == 2 || topo[p_topo_ID].type == 5) && i+h_c < clust_topo.size())
        {
            //cerr << "Tail-middle [" << i+h_c << "] == [" << p_topo_ID << "] type2|5==" << topo[ p_topo_ID ].type << endl;
            clust_topo[i+h_c] = topo[ p_topo_ID ]; // write tail-middle
        }
    }

    void add_lipid_last_bead(Atoms& topo, Atoms& clust_topo, size_t& i, size_t& h_c, int p_topo_ID)
    {
        if( (topo[p_topo_ID].type == 3 || topo[p_topo_ID].type == 6) && i+h_c < clust_topo.size() )
        {
            //cerr << "Tail-Last [" << i+h_c << "] == [" << p_topo_ID << "] type3|6==" << topo[ p_topo_ID ].type << endl;
            clust_topo[i+h_c] = topo[ p_topo_ID ]; // write tail-end
        }
    }

    Tensor_xyz get_dir(Data& data, Tensor_xyz pos_a, Tensor_xyz pos_b)
    {
        Tensor_xyz periodic = Tensor_xyz(data.in.sim_box.xhi - data.in.sim_box.xlo, data.in.sim_box.yhi - data.in.sim_box.ylo, data.in.sim_box.zhi - data.in.sim_box.zlo);
        Tensor_xyz dir = (pos_a - pos_b);

        // fix for periodic box
        if(dir.x >  0.5*periodic.x) { dir.x -= periodic.x; }
        if(dir.x < -0.5*periodic.x) { dir.x += periodic.x; }
        if(dir.y >  0.5*periodic.y) { dir.y -= periodic.y; }
        if(dir.y < -0.5*periodic.y) { dir.y += periodic.y; }
        if(dir.z >  0.5*periodic.z) { dir.z -= periodic.z; }
        if(dir.z < -0.5*periodic.z) { dir.z += periodic.z; }

        dir.normalise();

        return dir;
    }

    void get_dirs(Data& data, Atoms& frame, vector<Tensor_xyz>& dirs)
    {
        for(size_t i=0; i<frame.size(); i+=4)
        {
            dirs[i/4] = get_dir(data, frame[i].pos, frame[i+3].pos);
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

    Histogram_Spherical analyze_dirs(Data& data, Atoms& topo, Trajectory& traj, int number_of_averaged_frames)
    {
        vector<Tensor_xyz> dirs(topo.size() / 4, Tensor_xyz(0.0, 0.0, 0.0));
        get_frame_averaged_dirs(data, dirs, topo, traj, number_of_averaged_frames);

        Histogram_Spherical h_sp(data.in.p_vec_int["Histo_spherical_settings"][0], data.in.p_vec_int["Histo_spherical_settings"][1]);

        for(Tensor_xyz a : dirs)
        {
            h_sp.add(a, Tensor_xyz(0.0,0.0,1.0), 0.0);
        }

        return h_sp;

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


    unordered_set<int> get_bilayer_lipid_head_set(Data& data, Atoms& topo)
    {
        unordered_set<int> lipid_head_set;
        Cell_List cell_list;

        int paired_lipid_head_ID = -1;
        cell_list.init( data.in.sim_box ); // allocate memory
        cell_list.add( topo ); // populate the cell list

        for(size_t i=0; i<topo.size(); i+=4) // loop over lipis
        {
            paired_lipid_head_ID = get_most_antiparallel_lip_head_ID(data, topo, cell_list, i);

            if(paired_lipid_head_ID > -1)
            {
                lipid_head_set.insert(i); // add lipids head, if i already in set, its ignored
                lipid_head_set.insert(paired_lipid_head_ID);
            }
        }

        cell_list.delete_all();

        return lipid_head_set;
    }

    int get_most_antiparallel_lip_head_ID(Data& data, Atoms& topo, Cell_List& cell_list, size_t i)
    {
        double angle_max = 150;
        double cutoff = 2.0;
        cutoff *= cutoff;

        int paired_lipid_head_ID = -1;
        Atom& lip_head=topo[i];
        Atom& lip_tail=topo[i+3];

        Tensor_xyz dir_a, dir_b;

        if(lip_head.type == 4 && lip_tail.type == 6) // select helper lipids
        {
            dir_a = get_dir(data, lip_head.pos, lip_tail.pos); // lip dir
            cell_list.set_neighbors_pbc( lip_tail.pos ); // particles close to tail

            for(size_t j=0; j<cell_list.neighbors.size(); ++j) // Loop over cells
            {
                for(size_t nei_ID : (*cell_list.neighbors[j]) ) // Loop over neigbor particles in those cells
                {
                    Atom& other_tail=topo[nei_ID];
                    if(nei_ID-3 >= 0 && other_tail.type == 6 && lip_tail.pos.distSQ_pbc(other_tail.pos, cell_list.pbc, cell_list.pbc_inv) < cutoff) // neighbor is last bead of tail
                    {
                        Atom& other_head=topo[nei_ID-3];
                        dir_b = get_dir(data, other_head.pos, other_tail.pos);

                        if(dir_a.angleToDegrees(dir_b) > angle_max) // find the neighbor with angle above limit, select only the most antiparallel lipid
                        {
                            paired_lipid_head_ID = nei_ID-3;
                        }
                    }
                }
            }
        }

        return paired_lipid_head_ID;
    }
};


#endif // LIPID_NANOPARTICLE_H
