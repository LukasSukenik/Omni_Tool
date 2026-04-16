#ifndef RDF_H
#define RDF_H

#include <string>

#include "histogram.h"
#include "atom.h"

using namespace std;




void print_rdf(string filename, Histogram& histo)
{
    fstream fs(filename, fstream::out);
    if(!fs.is_open())
    {
        cout << "Failed to open file " << filename << endl;;
    }
    else
    {
        double r_start;
        double r_stop;
        const double pi = 3.141592653589793;
        double vol;
        for(size_t i=0; i<histo.size(); ++i)
        {
            r_start = histo.min+ ( ((histo.max-histo.min)*i) / histo.size() );
            r_stop = histo.min+ ( ((histo.max-histo.min)*(i+1)) / histo.size() );
            vol = ((4.0/3.0)*pi) * ((pow(r_stop,3.0))-(pow(r_start,3.0)));
            fs << r_start  << " " << histo.get(i) / vol  << endl;
        }
    }
    fs.close();
}




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
    print_rdf(outfile, hist);
}



#endif // RDF_H
