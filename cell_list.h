#ifndef CELL_LIST_H
#define CELL_LIST_H

#include <vector>
#include <cmath>
#include <stack>

#include "atom.h"
#include "sim_box.h"

using namespace std;

class Tensor_xyz_integer
{
public:
    Tensor_xyz_integer(){}
    Tensor_xyz_integer(int x, int y, int z) : x(x), y(y), z(z) {}
    int x,y,z;
};

class Lattice_3D
{
public:
    vector<vector<vector<bool>>> occupied;
    vector<vector<vector<bool>>> visited;

    const vector<Tensor_xyz_integer> directions = { {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1} };

    Lattice_3D(Tensor_xyz_integer cell_count)
    {
        occupied.resize(cell_count.x, vector<vector<bool>>(cell_count.y, vector<bool>(cell_count.z, false)));
        visited.resize(cell_count.x, vector<vector<bool>>(cell_count.y, vector<bool>(cell_count.z, false)));
    }

    void set_Occupied(Tensor_xyz_integer i, bool value = true)
    {
        if(in_Bounds(i)) {
            occupied[i.x][i.y][i.z] = value;
        }
    }

    void reset_occupied()
    {
        set_all(occupied, false);
    }

    void reset_visited()
    {
        set_all(visited, false);
    }

    void add_particle(Tensor_xyz& pos, Tensor_xyz& box, int number_of_cells, double inv_cell_size)
    {
        Tensor_xyz_integer i;
        i.x = binning_fce(pos.x, box.x, number_of_cells, inv_cell_size); // index of cell, where particle COM hit
        i.y = binning_fce(pos.y, box.x, number_of_cells, inv_cell_size);
        i.z = binning_fce(pos.z, box.x, number_of_cells, inv_cell_size);
        set_Occupied(i);
    }

    void dfs(Tensor_xyz_integer start)
    {
        if(!in_Bounds(start))                   { cout << "Start cell is out of bounds" << endl; return; }
        if(occupied[start.x][start.y][start.z]) { cout << "Start cell is occupied" << endl; return; }
        reset_visited();

        stack<Tensor_xyz_integer> st;
        st.push(start);
        visited[start.x][start.y][start.z] = true;

        while(!st.empty()) {
            Tensor_xyz_integer current = st.top();
            st.pop();

            for(const Tensor_xyz_integer& dir : directions)
            {
                Tensor_xyz_integer next{current.x + dir.x, current.y + dir.y, current.z + dir.z};

                if( can_Visit(next) )
                {
                    visited[next.x][next.y][next.z] = true;
                    st.push(next);
                }
            }
        }
    }

private:
    bool in_Bounds(Tensor_xyz_integer i) const
    {
        return i.x >= 0 && i.x < occupied.size() && i.y >= 0 && i.y < occupied[0].size() && i.z >= 0 && i.z < occupied[0][0].size();
    }

    bool can_Visit(Tensor_xyz_integer i) const
    {
        return in_Bounds(i) && !occupied[i.x][i.y][i.z] && !visited[i.x][i.y][i.z];
    }

    //
    // values are distributed from -box_size/2 to box_size/2
    //
    inline int binning_fce(double value, double box_size, int number_of_cells, double inv_cell_size)
    {
        if(value < -0.5*box_size) return 0;
        if(value >= 0.5*box_size) return number_of_cells - 1;

        return (value + 0.5*box_size) * inv_cell_size; // value / cell_size, where cell_size = box/number_of_cells
    }

    void set_all(vector<vector<vector<bool>>>& lattice, bool value=false)
    {
        for(auto& plane : lattice)
        {
            for(auto& row : plane)
            {
                fill(row.begin(), row.end(), value);
            }
        }
    }
};

/**
 * @brief The Cell_List class
 * Usage example in Lipid::populate()
 *
 * cell_list.init( data ); // allocate memory
 * coll_cell_list.add( Atoms ); // populate the cell list
 *
 * // calculate particle distance
 * cell_list.set_neighbors( part.get_center_of_mass().pos ); // setup particle for distance calc
 * cell_list.neighbors; //  is now the neighbor list of particle part
 *
 * cell_list.delete_all();
 *
 */
class Cell_List : public vector<vector<vector< vector<int >* >>> // 3D spatial indices + list of index to coll_beads
{
public:
    Cell_List(){}

    Tensor_xyz_integer number_of_cells;
    Tensor_xyz cell_size; // for lipids at leat 4.0
    Tensor_xyz pbc;
    Tensor_xyz pbc_inv;

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

    double round_to(double a, double decimals)
    {
        return round( a*pow(10.0,decimals) ) / pow(10.0,decimals);
    }

    /**
     * @brief init - allocate cell list memory
     */
    void init(Simulation_Box& sim_box)
    {
        double decimal_places = 3.0; // XTC is usually rounded to 3 decimal places, meaning with lammps 4 decimal places a particle can be out of bound
        double tiny_num = pow(10.0, -5.0);
        xlo = round_to(sim_box.xlo, decimal_places) -tiny_num; // just rounding is not enough, need to enlarge the box slightly
        xhi = round_to(sim_box.xhi, decimal_places) +tiny_num;
        ylo = round_to(sim_box.ylo, decimal_places) -tiny_num;
        yhi = round_to(sim_box.yhi, decimal_places) +tiny_num;
        zlo = round_to(sim_box.zlo, decimal_places) -tiny_num;
        zhi = round_to(sim_box.zhi, decimal_places) +tiny_num;

        pbc.x = xhi - xlo;
        pbc.x = yhi - ylo;
        pbc.x = zhi - zlo;

        pbc_inv.x = 1.0 / pbc.x;
        pbc_inv.y = 1.0 / pbc.y;
        pbc_inv.z = 1.0 / pbc.z;

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
        //cerr << "Cell_List::init_cell_list" << endl;
        //show();
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

    void add(Atoms& a, int offset=0)
    {
        if(!a.empty())
        {
            for(size_t i=0; i<a.size(); ++i)
            {
                /*if( get_spatial_X(a[i].pos) <0 || get_spatial_Y(a[i].pos)<0 || get_spatial_Z(a[i].pos)<0 ) // testing out-of-bound spatial indexes caused by XTC low accuracy
                {
                    show();
                    cerr << "Size: " << this->size() << ", " << this->at(0).size() << ", " << this->at(0).at(0).size() << endl;
                    cerr << a[i].pos << endl;
                    cerr << "Cell_List::add [" << get_spatial_X(a[i].pos) << ", " << get_spatial_Y(a[i].pos) << ", " << get_spatial_Z(a[i].pos) << "] = " << i << endl;
                }*/
                this->at(get_spatial_X(a[i].pos)).at(get_spatial_Y(a[i].pos)).at(get_spatial_Z(a[i].pos))->push_back(offset + i); // does boundary checking
                //(*this)[get_spatial_X(a[i].pos)][get_spatial_Y(a[i].pos)][get_spatial_Z(a[i].pos)]->push_back(offset + i);

            }
        }
    }

    void set_neighbors_pbc(Tensor_xyz pos) // No PBC version
    {
        //cerr << "Cell_List::set_neighbors start" << endl;
        int x = get_spatial_X(pos);
        int y = get_spatial_Y(pos);
        int z = get_spatial_Z(pos);
        int ii, jj, kk;
        for(int i=0; i<3; ++i) // keep 0,1,2 because of the neighbors index
        {
            for(int j=0; j<3; ++j)
            {
                for(int k=0; k<3; ++k)
                {
                    ii = (i-1 + x + number_of_cells.x) % number_of_cells.x;
                    jj = (j-1 + y + number_of_cells.y) % number_of_cells.y;
                    kk = (k-1 + z + number_of_cells.z) % number_of_cells.z;

                    neighbors[(i*9) +(j*3) + k] = (*this)[ii][jj][kk];
                }
            }
        }
        //cerr << "Cell_List::set_neighbors end" << endl;
        //show_neighbors();
    }

    void set_neighbors(Tensor_xyz pos) // No PBC version
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

private:
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
};

#endif // CELL_LIST_H
