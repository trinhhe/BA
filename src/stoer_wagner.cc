#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/stoer_wagner_min_cut.hpp>

#include "benchmark.h"
#include "bitmap.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "timer.h"

using namespace std;

int Stoer_Wagner(const WGraph &g)
{
    Timer t;
    typedef boost::property<boost::edge_weight_t, int> EdgeWeightProp;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProp> Graph;
    Graph conn(g.num_nodes());
    for (auto i : g.vertices())
        for(auto j : g.out_neigh(i))
        {
            if(i < j.v)
                boost::add_edge(i, j.v, j.w, conn);
        }
    
    auto weights = get(boost::edge_weight, conn);
    t.Start();
    int res = boost::stoer_wagner_min_cut(conn, weights);
    t.Stop();
    cout << g.num_nodes() << ", " << g.num_edges() << ", "<<  t.Seconds() << ", " << res << "\n";
    return res;
}

void DummyPrint(const WGraph &g, size_t min_cut_value)
{
    cout << "min cut value: " << min_cut_value << endl;
}

bool DummyVerifier(const WGraph &g, size_t test_min)
{
    return true;
}

int main(int argc, char *argv[])
{
    CLApp cli(argc, argv, "stoer_wagner");
    if (!cli.ParseArgs())
        return -1;
    WeightedBuilder b(cli);
    WGraph g = b.MakeGraph();
    BenchmarkKernel(cli, g, Stoer_Wagner,DummyPrint, DummyVerifier);
    return 0;
}