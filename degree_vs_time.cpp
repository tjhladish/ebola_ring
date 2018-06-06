#include <Percolation_Sim.h>

double mean_deg(const vector<Node*> nodes) {
    double mean = 0.0;
    for (Node* n: nodes) mean += n->deg();
    return mean / nodes.size();
}

int main() {
    // Construct Network
    Network net("name", Network::Undirected);
    //net.small_world(10000, 16, 0.5);
    net.populate(10000);
    //net.erdos_renyi(16);
    //net.rand_connect_powerlaw(0.8, 44.0);
    net.rand_connect_exponential(0.0647);

    for (int i = 0; i < 100; i++){
        // Choose and run simulation
        Percolation_Sim sim(&net);
        sim.set_transmissibility(0.075);
        sim.rand_infect(10);
        int time = 0;
        while (sim.count_infected() > 0) {
            cerr << i << " " << ++time << " " <<  sim.count_infected() << " " << mean_deg(sim.get_infected_nodes()) << endl;
            sim.step_simulation();
        }
        //sim.run_simulation();
        //cout << sim.epidemic_size() << endl;
        sim.reset();
    }
    return 0;
}
