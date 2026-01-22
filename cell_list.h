#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <vector>

#include "atom.h"
#include "data.h"

using namespace std;

class Tensor_xyz_integer
{
public:
    Tensor_xyz_integer(){}
    int x,y,z;
};

class Cell_List : public vector<vector<vector< vector<int >* >>> // 3D spatial indices + list of index to coll_beads
{
public:
    Cell_List(){}

    Tensor_xyz_integer number_of_cells;
    Tensor_xyz cell_size; // for lipids at leat 4.0

    double xlo = 0.0; // for spatial indices
    double xhi = 0.0;
    double ylo = 0.0;
    double yhi = 0.0;
    double zlo = 0.0;
    double zhi = 0.0;

    vector< vector<int>* > neighbors;
    vector<int>* empty_list;

    void show()
    {
        cerr << "xlo:" << xlo << " xhi:" << xhi << endl;
        cerr << "ylo:" << ylo << " yhi:" << yhi << endl;
        cerr << "zlo:" << zlo << " zhi:" << zhi << endl;
        cerr << "number of cells: " << number_of_cells.x << ", " << number_of_cells.y << ", " << number_of_cells.z << endl;
        cerr << "cell_size: " << cell_size.x << ", " << cell_size.y << ", " << cell_size.x << endl;
    }


    void init(Data& data)
    {
        xlo = data.in.sim_box.xlo;
        xhi = data.in.sim_box.xhi;
        ylo = data.in.sim_box.ylo;
        yhi = data.in.sim_box.yhi;
        zlo = data.in.sim_box.zlo;
        zhi = data.in.sim_box.zhi;

        number_of_cells.x = floor( (xhi - xlo)/4.0 );
        number_of_cells.y = floor( (yhi - ylo)/4.0 );
        number_of_cells.z = floor( (zhi - zlo)/4.0 );

        cell_size.x = (xhi - xlo) / number_of_cells.x;
        cell_size.y = (xhi - xlo) / number_of_cells.x;
        cell_size.z = (xhi - xlo) / number_of_cells.x;

        empty_list = new vector<int>();

        neighbors.resize(27);
        neighbors.assign(27, empty_list);

        (*this).resize( number_of_cells.x );

        for(int i=0; i<number_of_cells.x; ++i)
        {
            (*this)[i].resize(number_of_cells.y);

            for(int j=0; j<number_of_cells.y; ++j)
            {
                (*this)[i][j].resize(number_of_cells.z);
                for(int k=0; k<number_of_cells.z; ++k)
                {
                    (*this)[i][j][k] = new vector<int >();
                    (*this)[i][j][k]->reserve(100);
                }
            }
        }
        cerr << "Cell_List::init_cell_list" << endl;
        show();
    }

    void delete_all()
    {
        int number_of_cells = this->size();
        for(int i=0; i<number_of_cells; ++i)
        {
            for(int j=0; j<number_of_cells; ++j)
            {
                for(int k=0; k<number_of_cells; ++k)
                {
                    delete (*this)[i][j][k];
                }
            }
        }
        delete empty_list;
    }

    int get_spatial_X(Tensor_xyz pos)
    {
        return floor(floor(pos.x - xlo) / cell_size.x);
    }

    int get_spatial_Y(Tensor_xyz pos)
    {
        return floor(floor(pos.y - ylo) / cell_size.y);
    }

    int get_spatial_Z(Tensor_xyz pos)
    {
        return floor(floor(pos.z - zlo) / cell_size.z);
    }

    void add(Atoms& coll, int offset=0)
    {
        if(!coll.empty())
        {
            for(unsigned int i=0; i<coll.size(); ++i)
            {
                (*this)[get_spatial_X(coll[i].pos)]
                       [get_spatial_Y(coll[i].pos)]
                       [get_spatial_Z(coll[i].pos)]->push_back(offset + i);
                //cerr << "Cell_List::add_lipid_cell_list [" << get_spatial_index(part[0].pos.x) << ", " << get_spatial_index(part[0].pos.y) << ", " << get_spatial_index(part[0].pos.z) << "] = " << beads.size() + i << endl;
            }
        }
    }

    void set_neighbors(Tensor_xyz pos)
    {
        //cerr << "Cell_List::set_neighbors start" << endl;
        int x = get_spatial_X(pos);
        int y = get_spatial_Y(pos);
        int z = get_spatial_Z(pos);
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                for(int k=0; k<3; ++k)
                {
                    if(x+i-1 < 0 || y+j-1 < 0 || z+k-1 < 0 || x+i-1 >= number_of_cells.x || y+j-1 >= number_of_cells.y || z+k-1 >= number_of_cells.z )
                    {
                        neighbors[(i*9) +(j*3) + k] = empty_list;
                    }
                    else
                    {
                        neighbors[(i*9) +(j*3) + k] = (*this)[x+i-1][y+j-1][z+k-1];
                    }
                }
            }
        }
        //cerr << "Cell_List::set_neighbors end" << endl;
        //show_neighbors();
    }

    void show_neighbors()
    {
        int sum=0;
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                for(int k=0; k<3; ++k)
                {
                    cerr << "Neighbor:" << i << ", " << j << ", " << k << " : " << neighbors[(i*9) +(j*3) + k]->size() << endl;
                    sum += neighbors[(i*9) +(j*3) + k]->size();
                }
            }
        }
        cerr << "Neighbor sum: " << sum << endl;
    }

    bool validate_neighbors(Atom test, Atoms& beads)
    {
        double max_dist_SQ = 4.0 * ( cell_size.x*cell_size.x + cell_size.y*cell_size.y + cell_size.z*cell_size.z );

        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                for(int k=0; k<3; ++k)
                {
                    for(int index : (*neighbors[(i*9) +(j*3) + k]) )
                    {
                        if(test.distSQ(beads[index]) > max_dist_SQ)
                        {
                            cerr << "validate_neighbors::Maximal distance exceeded: " << "[" << i << ", " << j << ", " << k << "] : " << sqrt(test.distSQ(beads[index])) << " > " << sqrt(max_dist_SQ) << endl;
                            cerr << "    test[" << test.pos.x << ", " << test.pos.y << ", " << test.pos.z << "] vs beads[" << index << "]:[" << beads[index].pos.x << ", " << beads[index].pos.y << ", " << beads[index].pos.z << "]" << endl;
                            exit(1);
                        }
                    }
                }
            }
        }
        return true;
    }
};

#endif // CELL_LIST_H
