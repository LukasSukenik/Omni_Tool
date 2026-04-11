#ifndef RDF_H
#define RDF_H

#include <string>

#include "histogram.h"
#include "atom.h"

using namespace std;


void rdf(Atoms& topo, double min, double max, int h_size, string outfile)
{
    Histogram hist(h_size, min, max);
    Atom com = topo.get_center_of_mass();
    double r=0.0;
    for(Atom& a : topo)
    {
        r = a.pos.dist(com.pos);
        if(r > max)
        {
            cerr << "Vesicle sigma radius exeeds bounds of histogram set at " << max << ", r calculated at " << r << endl;
            exit(1);
        }
        hist.add(r);
    }
    hist.print(outfile);
}

#endif // RDF_H
