// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <string.h>
#include <assert.h>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "spantree.h"

/*
GAP Benchmark Suite
Kernel: MinCut
Author: Henry Trinh

Will output mincut

Requires input graph:
  - to be undirected
  - no duplicate edges 
  - neighborhoods are sorted by vertex identifiers

Other than symmetrizing, the rest of the requirements are done by SquishCSR
during graph building.

*/

using namespace std;
typedef EdgePair<NodeID, WNode> WEdge;

NodeID find(int *parent, NodeID i)
{
    NodeID root = i;
    while (root != parent[root])
        root = parent[root];
    // path compression
    while (i != root)
    {
        int newi = parent[i];
        parent[i] = root;
        i = newi;
    }
    return root;
}

void Union(int *parent, int *size, NodeID x, NodeID y)
{
    //union by size
    if (size[x] < size[y])
    {
        parent[x] = y;
        size[y] += size[x];
    }
    else
    {
        parent[y] = x;
        size[x] += size[y];
    }
}

pvector<WEdge> Kruskal(const WGraph &g)
{
    pvector<WEdge> WEdgelist(g.num_edges());
    int j = 0;
    // cout << "g.out_neigh(0)[1]: " << *(g.out_neigh(0).begin() + 3) << endl;
    for (NodeID u : g.vertices())
    {
        for (WNode wn : g.out_neigh(u))
        {
            //since g.vertices and g.out_neigh give us sorted NodeID, with u < wn.v we don't count edges twice
            if (u < wn.v)
            {
                WEdgelist[j] = WEdge(u, wn);
                j++;
            }
        }
    }
    // less function for sorting WEdgelist
    auto lessWEdge = [](WEdge const &l, WEdge const &r) {
        return l.v.w < r.v.w;
    };

    sort(WEdgelist.begin(), WEdgelist.end(), lessWEdge);
    //will contain edges from span tree
    pvector<WEdge> tree_WEdgelist(g.num_nodes() - 1);
    //array to help resolving cycles
    int *parent = new int[g.num_nodes()];
    //array to help for union by size
    int *size = new int[g.num_nodes()];
    for (int i = 0; i < g.num_nodes(); i++)
    {
        parent[i] = i;
        size[i] = 1;
    }

    int tree_edges = 0;
    // current edge from graph
    int i = 0;
    while (tree_edges < g.num_nodes() - 1)
    {
        WEdge next_edge = WEdgelist[i];
        NodeID x = find(parent, next_edge.u);
        NodeID y = find(parent, next_edge.v.v);
        if (x != y)
        {
            tree_WEdgelist[tree_edges] = next_edge;
            Union(parent, size, x, y);
            tree_edges++;
        }
        i++;
    }
    // for (auto i : tree_WEdgelist)
    //     cout << "Edges " << i.u << " " << i.v.v << " " << i.v.w << endl;

    delete[] parent;
    delete[] size;

    // cout << "#treeEdges: " << tree_edges << " " << endl;
    return tree_WEdgelist;
}

size_t MinCut(const WGraph &g)
{
    pvector<WEdge> tree_edges = Kruskal(g);

    // for (auto i : vertex_degree)
    // {
    //     cout << i << endl;
    // }

    auto tree_graph = WeightedBuilder::Load_CSR_From_Edgelist(tree_edges, true);
    SpanTree T(tree_graph);
    T.print();
    return 0;
}

void PrintMinCutValue(const WGraph &g, size_t min_cut_value)
{
    cout << "min cut value: " << min_cut_value << endl;
}

// Compares with simple serial implementation (Stoer-Wagner), not really efficient with n * n matrix, only for small samples for the sake of debugging
bool MINCUTVerifier(const WGraph &g, size_t test_min)
{
    int n = g.num_nodes();
    pvector<int> edges(n * n, 0);
    //
    int v[n], dis[n], vis[n];
    //construct edges in matrix form
    for (NodeID u : g.vertices())
    {
        assert(g.out_degree(u) > 0);
        // if (g.out_degree(u) == 0)
        // {
        //     cout << u << " has no neighbors" << endl;
        //     return false;
        // }
        for (WNode wn : g.out_neigh(u))
        {
            edges[u * n + wn.v] = wn.w;
            edges[wn.v * n + u] = wn.w;
        }
    }

    size_t mincut = numeric_limits<size_t>::max();
    for (int i = 0; i < n; i++)
    {
        v[i] = i;
        vis[i] = 0;
    }
    int n_ = n - 1;
    //loop until only 2 nodes left,
    while (n_ > 0)
    {
        //last and second last node, we start with node 1 and 0
        int p = 1, prev = 0;
        for (int i = 0; i <= n_; i++)
        {
            dis[v[i]] = edges[v[0] * n + v[i]];
            if (dis[v[i]] > dis[v[p]])
                p = i;
        }
        vis[v[0]] = n_;
        for (int i = 1; i <= n_; i++)
        {
            if (i == n_)
            {
                mincut = min<size_t>(mincut, dis[v[p]]);
                //contract edges from p and prev
                for (int j = 0; j <= n_; j++)
                {
                    edges[v[prev] * n + v[j]] += edges[v[p] * n + v[j]];
                    edges[v[j] * n + v[prev]] = edges[v[prev] * n + v[j]];
                }
                v[p] = v[n_--];
                break;
            }
            vis[v[p]] = n_;
            prev = p;
            p = -1;
            for (int j = 1; j <= n_; j++)
            {
                if (vis[v[j]] != n_)
                {
                    //find most tightly connected vertex which hasn't been visited yet
                    dis[v[j]] += edges[v[prev] * n + v[j]];
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
    // g.PrintTopology();
    if (g.directed())
    {
        cout << "Input graph is directed but we only consider undirected graphs" << endl;
        return -2;
    }
    BenchmarkKernel(cli, g, MinCut, PrintMinCutValue, MINCUTVerifier);
}
