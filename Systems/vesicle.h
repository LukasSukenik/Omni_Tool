#ifndef VESICLE_H
#define VESICLE_H

#include "system_base.h"
#include "atom.h"
#include "xtcanalysis.h"


#include "sphere.h"
#include "lipid.h"

class Vesicle : public System_Base, public Sphere
{
public:
    inline static const string keyword = "Vesicle";
    const string name = "Vesicle";

    Vesicle() : System_Base("Vesicle") {}

    void execute(Data& data)
    {
        int sys_id = data.id_map[data.in.id];
        Atoms& ves = data.coll_beads[sys_id];

        if(data.in.system_function.compare("rdf") == 0)
        {
            radial_distribution_histogram(ves, data.in.system_var_a, data.in.system_var_b);
            cerr << "Generated radial distribution function" << endl;
        }

        if(data.in.system_function.compare("Make_2_Vesicle_System") == 0)
        {
            data.coll_beads.push_back(ves);
            Atoms& ves = data.coll_beads[ sys_id ]; // push_back can allocate memory, which will brick the reference ves
            Atoms& ves2 = data.coll_beads.back();

            ves2.set_mol_tag(2);
            while( ves.is_overlap(ves2) )
            {
                ves.move(Atom(1.0, 0.0, 0.0));
                ves2.move(Atom(-1.0, 0.0, 0.0));
            }

            cerr << "Generated 2 vesicle system" << endl;
        }
    }

    void analyze_leaflet_lipid_count(Atoms& ves)
    {
        Lipids vesicle;
        for(int i=0; i<ves.size(); i+=4)
        {
            vesicle.push_back( Lipid(ves[i], ves[i+1], ves[i+2], ves[i+3]) );
        }

        Trajectory traj;
        traj.load("file.xtc");

        int c_up=0;
        int c_low=0;
        get_leaflet_count(vesicle, c_up, c_low);
        cout << "0 " << c_up << " " << c_up << " " << c_low << " " << c_low << endl;

        int count_upper=0;
        int count_lower=0;
        for(unsigned int i=0; i<traj.size(); ++i)
        {
            vector<Tensor_xyz>& frame = traj[i];
            if( vesicle.update_positions(frame) )
            {
                count_upper=0;
                count_lower=0;
                get_leaflet_count(vesicle, count_upper, count_lower);
                cout << i+1 << " " << count_upper << " " << c_up << " " << count_lower << " " << c_low << endl;
            }
        }
    }

    void radial_distribution_histogram(Atoms& at_col, double r_max, int hist_size)
    {
        double r=0.0;
        vector<int> histogram(hist_size,0);
        Atom com = at_col.get_center_of_mass();

        for(Atom& a : at_col)
        {
            r = a.pos.dist(com.pos);
            if(r > r_max)
            {
                cerr << "Vesicle sigma radius exeeds bounds of histogram set at " << r_max << ", r calculated at " << r << endl;
                exit(1);
            }
            histogram[(int) hist_size * (r/r_max)]++;
        }

        for(unsigned int i=0; i<histogram.size(); ++i)
        {
            cout << r_max * i / histogram.size() << " " << histogram[i] << endl;
        }
    }

    void radial_distribution_histogram(vector<Tensor_xyz>& frame)
    {
        double r_max=20;
        double r=0.0;
        vector<int> histogram(1000,0);
        Tensor_xyz com(0.0,0.0,0.0);

        for(Tensor_xyz& pos : frame)
        {
            com += pos;
        }
        com /= frame.size();

        for(Tensor_xyz& pos : frame)
        {
            r = pos.dist(com);
            histogram[(int) 1000 * (r/r_max)]++;
        }

        for(unsigned int i=0; i<histogram.size(); ++i)
        {
            cout << r_max * i / histogram.size() << " " << histogram[i] << endl;
        }
    }

    void get_leaflet_count(Lipids& vesicle, int& count_upper, int& count_lower)
    {
        Atom dir;
        Atom head_dir;

        for(Lipid& lip : vesicle)
        {
            dir = lip.get_direction();
            head_dir = lip.part[0];
            head_dir.normalise();

            if( dir.dot(head_dir) > 0.0)
            {
                count_upper++;
            }
            else
            {
                count_lower++;
            }
        }
    }


    void generate( Data& data )
    {
        test_input(data);

        double outer_leaflet_radius = data.in.radius+4.0;
        double inner_leaflet_area = pow(data.in.radius, 2);
        double outer_leaflet_area = pow(outer_leaflet_radius, 2);
        double sphere_surface_relative_increase = outer_leaflet_area / (inner_leaflet_area);


        int lower_leaflet_lipid_count = data.in.num_lipids * 1.0 / (1.0 + sphere_surface_relative_increase);
        int lower_leaflet_receptor_count = data.in.num_rec * 1.0/ (1.0 + sphere_surface_relative_increase);
        int upper_leaflet_lipid_count = data.in.num_lipids * sphere_surface_relative_increase / (1.0 + sphere_surface_relative_increase);
        int upper_leaflet_receptor_count = data.in.num_rec * sphere_surface_relative_increase / (1.0 + sphere_surface_relative_increase);

        vector<Lipid> lower_leaflet = gen_lower_leaflet(lower_leaflet_lipid_count, data.in.mol_tag, data.in.radius);
        convert_receptors(lower_leaflet, lower_leaflet_receptor_count);

        vector<Lipid> upper_leaflet = gen_upper_leaflet(upper_leaflet_lipid_count, data.in.mol_tag, outer_leaflet_radius, lower_leaflet.size()*4);
        convert_receptors(upper_leaflet, upper_leaflet_receptor_count);

        for(Lipid& lip : lower_leaflet)
        {
            beads.insert(beads.end(), lip.part.begin(), lip.part.end());
            bonds.insert(bonds.end(), lip.bond.begin(), lip.bond.end());
        }

        for(Lipid& lip : upper_leaflet)
        {
            beads.insert(beads.end(), lip.part.begin(), lip.part.end());
            bonds.insert(bonds.end(), lip.bond.begin(), lip.bond.end());
        }

        if(lower_leaflet.empty() || upper_leaflet.empty())
        {
            cerr << "Vesicle::generate leaflet empty" << endl;
            exit(1);
        }
    }

    void test_input(Data& data)
    {
        if(data.in.num_rec > data.in.num_lipids)
        {
            cerr << "Number of receptors " << data.in.num_rec << " cannot be larger than number of lipids " << data.in.num_lipids << endl;
            exit(1);
        }
    }

    vector<Lipid> gen_lower_leaflet(int num_lipids, int mol_tag, double radius)
    {
        vector<Lipid> lipids;

        Atoms hd = fibonacci_sphere_radius(num_lipids, mol_tag, radius     );
        Atoms t1 = fibonacci_sphere_radius(num_lipids, mol_tag, radius +1.0);
        Atoms t2 = fibonacci_sphere_radius(num_lipids, mol_tag, radius +2.0);
        Atoms t3 = fibonacci_sphere_radius(num_lipids, mol_tag, radius +3.0);

        for(unsigned int i=0; i<hd.size(); ++i)
        {
            hd[i].N = 4*i+1;
            t1[i].N = 4*i+2;
            t2[i].N = 4*i+3;
            t3[i].N = 4*i+4;

            lipids.push_back( Lipid(hd[i], t1[i], t2[i], t3[i]) );
            lipids[i].set_bead_type(Lipid::Leaflet::upper);
        }

        return lipids;
    }

    vector<Lipid> gen_upper_leaflet(int num_lipids, int mol_tag, double radius, int offset)
    {
        vector<Lipid> lipids;

        Atoms hd = fibonacci_sphere_radius(num_lipids, mol_tag, radius +3.0);
        Atoms t1 = fibonacci_sphere_radius(num_lipids, mol_tag, radius +2.0);
        Atoms t2 = fibonacci_sphere_radius(num_lipids, mol_tag, radius +1.0);
        Atoms t3 = fibonacci_sphere_radius(num_lipids, mol_tag, radius     );

        for(unsigned int i=0; i<hd.size(); ++i)
        {
            hd[i].N = 4*i+1 + offset;
            t1[i].N = 4*i+2 + offset;
            t2[i].N = 4*i+3 + offset;
            t3[i].N = 4*i+4 + offset;

            lipids.push_back( Lipid(hd[i], t1[i], t2[i], t3[i]) );
            lipids[i].set_bead_type(Lipid::Leaflet::upper);
        }

        return lipids;
    }

    Atoms fibonacci_sphere_radius(int num_lipids, int mol_tag, double radius)
    {
        Atoms a = fibonacci_sphere( num_lipids, 0, mol_tag);
        a.scale(radius);
        return a;
    }

    void convert_receptors(vector<Lipid>& lipids, int receptor_count)
    {
        int lipid_count = lipids.size();
        int count=0;
        int random=0;
        while(count < receptor_count)
        {
            random = (int)(ran() * lipid_count); // select random lipid
            if(random < lipids.size())
            {
                if(lipids[random].part[0].type == Lipid::head_upper_leaf)
                {
                    lipids[random].part[0].type = Lipid::receptor;
                    count++;
                }
                if(lipids[random].part[0].type == Lipid::head_lower_leaf)
                {
                    lipids[random].part[0].type = Lipid::receptor;
                    count++;
                }
            }
        }
    }
};

#endif // VESICLE_H
