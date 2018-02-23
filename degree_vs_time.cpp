#include <Percolation_Sim.h>

double mean_deg(const vector<Node*> nodes) {
    double mean = 0.0;
    for (Node* n: nodes) mean += n->deg();
    return mean / nodes.size();
}

int main() {
    // Construct Network
    Network net("name", Network::Undirected);
    net.populate(10000);
    net.erdos_renyi(20);

    for (int i = 0; i < 1; i++){
        // Choose and run simulation
        Percolation_Sim sim(&net);
        sim.set_transmissibility(0.075);
        sim.rand_infect(10);
        cerr << sim.count_infected() << " " << mean_deg(sim.get_infected_nodes()) << endl;
        while (sim.count_infected() > 0) {
            sim.step_simulation();
            cerr << sim.count_infected() << " " << mean_deg(sim.get_infected_nodes()) << endl;
        }
        //sim.run_simulation();
        cout << sim.epidemic_size() << endl;
        sim.reset();
    }
    return 0;
}
