#ifndef SURFACE_H
#define SURFACE_H

#include "atom.h"

/**
 * @brief The Surface class
 *
 * Test if
 *
 */
class Surface
{
public:
    Surface(array<Atom,12>& edge, vector< Atom>& beads) : edge(edge) {}

    array<Atom,12>& edge;

    virtual bool operator() (vector<Atom>& container, int size, int i, int j, int a, int b, int c, Atom& push, vector<int>& types)
    {
        bool same = false;

        for(unsigned int q=0; q< container.size(); q++) {
            if(container[q].isAproxSame(push)) {
                same = true;
            }
        }

        return !same;
    }
};

class PentaSurface : public Surface
{
public:
    using Surface::edge;

    PentaSurface(array<Atom,12>& edge, vector< Atom>& beads) : Surface(edge, beads) {}

    bool operator() (vector<Atom>& container, int size, int i, int j, int a, int b, int c, Atom& push, vector<int>& types) override
    {
        bool smooth = true;
        bool same = false;
        double margin = -0.0001;

        Atom center(edge[a].pos.x + edge[b].pos.x + edge[c].pos.x, edge[a].pos.y + edge[b].pos.y + edge[c].pos.y, edge[a].pos.z + edge[b].pos.z + edge[c].pos.z);
        Atom ab = edge[a] + edge[b];
        Atom ac = edge[a] + edge[c];
        Atom bc = edge[b] + edge[c];

        center *= 1.0/3.0;
        ab *= 0.5;
        ac *= 0.5;
        bc *= 0.5;

        push.type = 12;

        if(size == 13) { // smooth vs keyed
            if(i==6 && j==0)  // midpoints, always split == move to vertice A (line on certices A and B)
                push = ( (edge[b] - edge[a]) * (-0.5/(size-1)) ) + push;
            if(i==6 && j==6)
                push = ( (edge[c] - edge[a]) * (-0.5/(size-1)) ) + push;
            if(i==12 && j==6)
                push = ( (edge[b] - edge[c]) * (-0.5/(size-1)) ) + push;

            if(smooth) {
                if(i==7 && j==2)
                    push = ( (edge[b] - edge[a]) * (-0.5/(size-1)) ) + push;
                if(i==7 && j==5)
                    push = ( (edge[c] - edge[a]) * (-0.5/(size-1)) ) + push;
                if(i==10 && j==5)
                    push = ( (edge[b] - edge[c]) * (-0.5/(size-1)) ) + push;
            }

            if(i==8 && j==4) { // center between 3 pantamers -> lets have a hole
                return false;
            }
        }

        //
        // SET TYPES FOR PENTAMERS
        //
        if( ( (edge[a]-push).pos.size() - margin < (edge[b]-push).pos.size() ) && ( (edge[a]-push).pos.size() - margin < (edge[c]-push).pos.size() ) )
            push.type = a;
        if( ( (edge[b]-push).pos.size() - margin < (edge[a]-push).pos.size() ) && ( (edge[b]-push).pos.size() - margin < (edge[c]-push).pos.size() ) )
            push.type = b;
        if( ( (edge[c]-push).pos.size() - margin < (edge[a]-push).pos.size() ) && ( (edge[c]-push).pos.size() - margin < (edge[b]-push).pos.size() ) )
            push.type = c;


        for(unsigned int q=0; q< container.size(); q++) {
            if(container[q].isAproxSame(push)) {
                same = true;
            }
        }

        return !same;
    }
};

#endif // SURFACE_H

