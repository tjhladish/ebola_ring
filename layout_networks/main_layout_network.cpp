#include <iostream>
#include "ForceLayout.h"
#include <random>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <assert.h>

using namespace std;

int main(int argc, char** argv) {
    default_random_engine rng;
    const int RES = 10000;
    uniform_real_distribution<double> runif(0.0,1.0);
    vector<Particle*> nodelist;

    const int W = 0.75*RES;
    const int H = 0.75*RES;

    assert(argc == 2);
    
    string buffer;
    istringstream line;
    ifstream fh_in(argv[1]);
    //stringstream ss;
    //ofstream fh_out("w_and_s_locs_w_id.txt");
    map<int, Particle*> nodes_lookup;
    if (fh_in) {
        cerr << "reading network\n";
        while (getline(fh_in, buffer)) {
            line.clear();
            line.str(buffer);
            int id1, id2;
            if (line >> id1 >> id2) {
                Particle* n1;
                Particle* n2;
                if (nodes_lookup.count(id1) == 0) {
                    n1 = new Particle(id1, 0, 0);
                    nodes_lookup[id1] = n1;
                    nodelist.push_back(n1);
                } else {
                    n1 = nodes_lookup[id1];
                }
                if (nodes_lookup.count(id2) == 0) {
                    n2 = new Particle(id2, runif(rng), runif(rng));
                    nodes_lookup[id2] = n2;
                    nodelist.push_back(n2);
                } else {
                    n2 = nodes_lookup[id2];
                }
                n1->addLinkTo(n2);
            }
        }
    }
    fh_in.close();

    
    for (unsigned int i = 0; i < nodelist.size(); ++i) {
        Particle* n = nodelist[i];
        if (i == 0) {
            n->setPos(W/2.0, H/2.0);
        } else {
            const double x = (runif(rng)*W - W/2)/(n->totalDegree() + 1);
            const double y = (runif(rng)*H - H/2)/(n->totalDegree() + 1);
            n->setPos(x, y);
        }
    }

/*
    nodelist.push_back(n1);
    for (int i = 0; i < 10; ++i) {
        Particle* n2 = new Particle(runif(rng), runif(rng));
        n1->addLinkTo(n2);
        nodelist.push_back(n2);
    }
    for (auto n: nodelist) {
        cerr << n->x() << " " << n->y() << endl;
    }
    cerr << endl;
 */   
    cerr << "laying out network\n";
    ForceLayout layout;
                        //L  R  B  T
    layout.set_dimensions(0, 10000, 0, 10000);
    int iterations = 1000;
    for (int t = 0; t < iterations; ++t) {
        layout.doLayout(nodelist, 1);
    }
     
    for (auto n: nodelist) {
        cout << n->get_id() << " " << n->x() << " " << n->y() << endl;
    }
}
