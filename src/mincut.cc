// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <string.h>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h" /*  */
#include "graph.h"
#include "pvector.h"

/*
GAP Benchmark Suite
Kernel: Triangle Counting (TC)
Author: Scott Beamer

Will count the number of triangles (cliques of size 3)

Requires input graph:
  - to be undirected
  - no duplicate edges (or else will be counted as multiple triangles)
  - neighborhoods are sorted by vertex identifiers

Other than symmetrizing, the rest of the requirements are done by SquishCSR
during graph building.

This implementation reduces the search space by counting each triangle only
once. A naive implementation will count the same triangle six times because
each of the three vertices (u, v, w) will count it in both ways. To count
a triangle only once, this implementation only counts a triangle if u > v > w.
Once the remaining unexamined neighbors identifiers get too big, it can break
out of the loop, but this requires that the neighbors to be sorted.

Another optimization this implementation has is to relabel the vertices by
degree. This is beneficial if the average degree is high enough and if the
degree distribution is sufficiently non-uniform. To decide whether or not
to relabel the graph, we use the heuristic in WorthRelabelling.
*/

using namespace std;

class spanTree
{
private:
    pvector<NodeID> parent;
    pvector<NodeID> depth; //useful when walking up the tree to find LCA for path decomposition
    pvector<NodeID> size;  //size of subtree rooted with pvector[i]
    pvector<NodeID> chain;
};

size_t MinCut(const WGraph &g)
{

    return (size_t)0;
}

void PrintMinCutValue(const WGraph &g, size_t min_cut_value)
{
    cout << "min cut value: " << min_cut_value << endl;
}

// Compares with simple serial implementation (Stoer-Wagner)
bool MINCUTVerifier(const WGraph &g, size_t test_min)
{
    int n = g.num_nodes();
    // int edges[n][n];
    int edges[n][n], v[n], dis[n], vis[n];
    memset(edges, 0, sizeof edges);
    for (NodeID u : g.vertices())
    {
        for (WNode wn : g.out_neigh(u))
        {
            edges[u][wn.v] = wn.w;
            edges[wn.v][u] = wn.w;
        }
    }

    size_t mincut = numeric_limits<size_t>::max();
    for (int i = 0; i < n; i++)
    {
        v[i] = i;
        vis[i] = 0;
    }
    int n_ = n - 1;
    while (n_ > 1)
    {
        int p = 1, prev = 0;
        for (int i = 1; i < n; i++)
        {
            dis[v[i]] = edges[v[0]][v[i]];
            if (dis[v[i]] > dis[v[p]])
                p = i;
        }
        vis[v[0]] = n_;
        for (int i = 1; i < n; i++)
        {
            if (i == n_)
            {
                mincut = min<size_t>(mincut, dis[v[p]]);
                //contract edges from p and prev
                for (int j = 0; j < n; j++)
                {
                    edges[v[prev]][v[j]] += edges[v[p]][v[j]];
                    edges[v[j]][v[prev]] = edges[v[prev]][v[j]];
                }
                v[p] = v[n_--];
                break;
            }
            vis[v[p]] = n_;
            prev = p;
            p = -1;
            for (int j = 1; j < n; j++)
            {
                if (vis[v[j]] != n_)
                {
                    dis[v[j]] += edges[v[prev]][v[j]];
                    if (p == -1 || dis[v[p]] < dis[v[j]])
                        p = j;
                }
            }
        }
    }
    if (mincut != test_min)
        cout << mincut << " != " << test_min << endl;
    return mincut == test_min;
}

int main(int argc, char *argv[])
{
    CLApp cli(argc, argv, "mincut");
    if (!cli.ParseArgs())
        return -1;
    WeightedBuilder b(cli);
    WGraph g = b.MakeGraph();
    if (g.directed())
    {
        cout << "Input graph is directed but we only consider undirected graphs" << endl;
        return -2;
    }
    BenchmarkKernel(cli, g, MinCut, PrintMinCutValue, MINCUTVerifier);
    return 0;
}
