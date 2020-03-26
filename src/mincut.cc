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
Kernel: MinCut
Author: Henry Trinh

Will output mincut

Requires input graph:
  - to be undirected
  - no duplicate edges (or else will be counted as multiple triangles)
  - neighborhoods are sorted by vertex identifiers

Other than symmetrizing, the rest of the requirements are done by SquishCSR
during graph building.

*/

using namespace std;
typedef EdgePair<NodeID, WNode> WEdge;

template <class = NodeID>
// class MinSpanTree
// {
// private:
//     pvector<NodeID> parent;
//     pvector<NodeID> size; //size of subtree rooted with pvector[i]
//     pvector<NodeID> chain;
//     // some edgelist here also
// public:
// };

int find(int *parent, int i)
{
    int root = i;
    while (root != parent[root])
        root = parent[root];
    //path compression
    while (i != root)
    {
        int newi = parent[i];
        parent[i] = root;
        i = newi;
    }
    return root;
}

void Union(int *parent, int *size, int x, int y)
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
    for (NodeID u : g.vertices())
    {
        for (WNode wn : g.out_neigh(u))
        {
            //since g.vertices and g.out_neigh give us sorted NodeID, with u < wn.v we don't account edges twice
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
    //edges from span tree
    pvector<WEdge> tree_WEdgelist(g.num_nodes() - 1);
    int *parent = new int[g.num_nodes()];
    int *size = new int[g.num_nodes()];
    for (int i = 0; i < g.num_nodes(); i++)
    {
        parent[i] = i;
        size[i] = 1;
    }

    int tree_edges = 0;
    // current edge from graph
    int i = 0;
    // cout << WEdgelist[i].u << " " << WEdgelist[i].v.v << " " << WEdgelist[i].v.w << endl;
    while (tree_edges < g.num_nodes() - 1)
    {
        WEdge next_edge = WEdgelist[i];
        int x = find(parent, next_edge.u);
        int y = find(parent, next_edge.v.v);
        if (x != y)
        {
            tree_WEdgelist[tree_edges] = next_edge;
            Union(parent, size, x, y);
            tree_edges++;
        }
        i++;
    }
    for (auto i : tree_WEdgelist)
        cout << "Edges " << i.u << " " << i.v.v << " " << i.v.w << endl;

    delete[] parent;
    delete[] size;

    cout << "#treeEdges: " << tree_edges << " " << endl;
    return tree_WEdgelist;
}

size_t MinCut(const WGraph &g)
{
    auto tree = Kruskal(g);
    return 0;
}

void PrintMinCutValue(const WGraph &g, size_t min_cut_value)
{
    cout << "min cut value: " << min_cut_value << endl;
}

// Compares with simple serial implementation (Stoer-Wagner)
bool MINCUTVerifier(const WGraph &g, size_t test_min)
{
    int n = g.num_nodes();
    pvector<int> edges(n * n, 0);
    //
    int v[n], dis[n], vis[n];
    //construct edges in matrix form
    for (NodeID u : g.vertices())
    {
        if (g.out_degree(u) == 0)
        {
            cout << u << " has no neighbors" << endl;
            return false;
        }
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
