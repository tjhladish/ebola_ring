#include "Network.h"
#include "Ring_Metrics.h"
#include <filesystem>

using namespace std;

vector<vector<int> > read_cols_of_ints(string filename, char sep=' ') {
 //   cerr << "Loading " << filename << endl;
    ifstream myfile(filename.c_str());

    vector<vector<int> > M;
    if (myfile.is_open()) {
        string line;

        while ( getline(myfile,line) ) {
            vector<string> fields;
            split(line, sep, fields);

            vector<int>row(fields.size());
            for( unsigned int i=0; i < fields.size(); i++ ) {
                    row[i] = stoi(fields[i]);
            }
            M.push_back(row);
        }
    }
    return M;
}

void process_level_data(Network &net, const vector<vector<int>>& level_data, vector<set<const Node*, NodePtrComp> > &levels, map<const Node*, int> &level_of) {
//int maxlvl = 0;
    for (auto row: level_data) {
//cerr << row[0] << " -- " << row[1] << endl;
        Node* n = net.get_node_by_name(to_string(row[0]));
        const int lvl = row[1];
//if (lvl > maxlvl) maxlvl = lvl;
        if ((signed) levels.size() <= lvl) levels.resize(lvl+1, set<const Node*, NodePtrComp>());

        levels[lvl].insert(n);
        level_of[n] = lvl;
    }
//cerr << "max lvl: " << maxlvl << endl;
}

void assign_interview_state(Network &net, vector<vector<int>> state_data) {
    for (auto row: state_data) {
// cerr << row[0] << ", " << row[1] << endl;
        net.get_node_by_name(to_string(row[0]))->set_state(row[1]);
    }
}

void usage() {
    cerr << "\n\tUsage: ./calculate_metrics RVT <network.csv> <levels.csv>\n\n";
    cerr << "\t       ./calculate_metrics GINMIX <network.csv <levels.csv> <interview.csv>\n\n";
}

int main(int argc, char** argv) {
    if (!(argc == 4 and string(argv[1]) == "RVT") and !(argc == 5 and string(argv[1]) == "GINMIX")) {usage(); exit(100);};
    string data_origin = argv[1];
    string net_file = argv[2];
    string lev_file = argv[3];

    assert(filesystem::exists(net_file));
    assert(filesystem::exists(lev_file));
    Network net(Network::Undirected);
    net.read_edgelist(net_file);

    /*
    cerr << "size: " << net.size() << endl;
    DistanceMatrix dm = net.get_node(0)->min_path_map();
    for (auto pair_nd: dm) {
        cerr << pair_nd.first->get_id() << " " << pair_nd.second << endl;
    }*/

    vector<vector<int>> level_data = read_cols_of_ints(lev_file);
    vector<set<const Node*, NodePtrComp> > levels;
    map<const Node*, int> level_of;

    process_level_data(net, level_data, levels, level_of);

    if (data_origin == "RVT") {
        TrialRawMetrics trm;
        raw_trial_metrics(levels, trm);
        trm.dumper(cout);
    } else if (data_origin == "GINMIX") {
        InterviewRawMetrics irm;
        string ivw_file = argv[4];
        assert(filesystem::exists(ivw_file));
        vector<vector<int>> interview_data = read_cols_of_ints(ivw_file);
        assign_interview_state(net, interview_data);
        PairwiseDistanceMatrix pdm = net.calculate_distances_map();

        raw_interview_metrics(&net, levels, level_of, pdm, irm);
        irm.dumper(cout);
    }
}
