// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <vector>
#include <string.h>
#include <assert.h>
#include <cmath>
#include <utility>
#include <map>
#include <set>

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "timer.h"
#include "util.h"
#include "spantree.h"
#include "minimumpath.h"

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

packing improvement idea 1 & 5
*/

using namespace std;
typedef EdgePair<NodeID, WNode> WEdge;
Timer t, t1, t2;

void dfs(WGraph &g, int v, pvector<bool> &visited)
{
    // visited[v] = true;
    // for (auto i : g.out_neigh(v))
    // {
    //     if (!visited[i])
    //         dfs(g, i, visited);
    // }
    stack<int> st;
    st.push(v);
    while (!st.empty())
    {
        auto s = st.top();
        st.pop();
        if (!visited[s])
            visited[s] = true;
        for (auto i : g.out_neigh(s))
        {
            if (!visited[i])
                st.push(i);
        }
    }
}

bool isConnected(WGraph &g)
{
    pvector<bool> visited(g.num_nodes(), false);
    dfs(g, 0, visited);
    for (int i : visited)
        if (!i)
            return false;
    return true;
}

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

//if min = true, we do MinST else MaxST
pvector<WEdge> Kruskal(const pvector<WEdge> &WEdgelist, int n, bool min)
{

    //copy to preserve WEdgelist order
    pvector<WEdge> edges(WEdgelist.begin(), WEdgelist.end());

    // less function for sorting WEdgelist
    auto lessWEdge = [](WEdge const &l, WEdge const &r) {
        return l.v.w < r.v.w;
    };

    auto moreWEdge = [](WEdge const &l, WEdge const &r) {
        return l.v.w > r.v.w;
    };

    if (min)
        sort(edges.begin(), edges.end(), lessWEdge);
    else
        sort(edges.begin(), edges.end(), moreWEdge);

    // //will contain edges from span tree
    pvector<WEdge> tree_WEdgelist(n - 1);
    //array to help resolving cycles
    int *parent = new int[n];
    //array to help for union by size
    int *size = new int[n];
    for (int i = 0; i < n; i++)
    {
        parent[i] = i;
        size[i] = 1;
    }

    int tree_edges = 0;
    // current edge from graph
    int i = 0;
    while (tree_edges < n - 1)
    {
        WEdge next_edge = edges[i];
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

    delete[] parent;
    delete[] size;

    return tree_WEdgelist;
}

//outputs MST as edgelist with edge id (sorted)
pvector<int> *KruskalWithLoad(const pvector<WEdge> &H, vector<pair<int, int>> &turn_number, int n)
{
    // pvector<pair<double, int>> turn_number_copy(turn_number.begin(), turn_number.end());

    sort(turn_number.begin(), turn_number.end());
    // contains the edge ids for mst
    pvector<int> *tree_list = new pvector<int>(n - 1);
    //array to help resolving cycles
    int *parent = new int[n];
    //array to help for union by size
    int *size = new int[n];
    for (int i = 0; i < n; i++)
    {
        parent[i] = i;
        size[i] = 1;
    }

    int tree_edges = 0;
    // current edge from graph
    size_t j = 0;

    while (tree_edges < n - 1)
    {
        // assert(turn_number[j].second >= 0);
        pair<int, int> edge = turn_number[j];
        NodeID x = find(parent, H[(size_t)edge.second].u);
        NodeID y = find(parent, H[(size_t)edge.second].v.v);
        if (x != y)
        {
            (*tree_list)[tree_edges] = edge.second;
            Union(parent, size, x, y);
            tree_edges++;
        }
        j++;
    }
    delete[] parent;
    delete[] size;
    return tree_list;
}

pair<double, pvector<pvector<int> *>> PackingWeight(
    const pvector<WEdge> &H, pvector<pair<int, int>> &remaining_capacity,
    double eps, int n, int m, int M, pvector<double> &tree_weights)
{
    //pair contains the edge id with its turn_number
    vector<pair<int, int>> turn_number;
    turn_number.reserve(m);
    double packing_value = 0;
    double inc = (eps * eps) / (3.0 * log(M));
    pair<double, pvector<pvector<int> *>> res;
    res.second.reserve(m + n * log(n) * log(n));
    for (int i = 0; i < m; i++)
        turn_number.push_back(make_pair(0, i));

    bool stop = false;
    // int iteration = 0;
    while (true)
    {
        pvector<int> *T = KruskalWithLoad(H, turn_number, n);

        // packing_value += inc;
        res.second.push_back(T);
        //update turn_number of edges contained in T
        int smallest_remaining_capacity = numeric_limits<int>::max();
        for (auto i : *T)
        {
            if (remaining_capacity[i].second < smallest_remaining_capacity)
                smallest_remaining_capacity = remaining_capacity[i].second;
        }
        for (auto i : *T)
        {
            remaining_capacity[i].second -= smallest_remaining_capacity;
            if (remaining_capacity[i].second == 0)
            {
                turn_number[i].first++;
                remaining_capacity[i].second = remaining_capacity[i].first;
                if (turn_number[i].first >= (3.0 * log(M)) / (eps * eps) - 1)
                    stop = true;
            }
        }
        tree_weights.push_back(smallest_remaining_capacity * inc);
        packing_value += smallest_remaining_capacity * inc;
        // iteration++;
        if (stop)
        {
            // cout << "packing iterations: " << iteration << endl;
            res.first = packing_value;
            return res;
        }
    }
}

// double binomial(int trials, int weight_cap, double p, default_random_engine gen)
// {
//     const double max_allowed_deviation = 1e-10;
//     if (p > 1 - max_allowed_deviation)
//         return min(trials, weight_cap);
//     uniform_real_distribution<double> distribution(0.0, 1.0);
//     double u = distribution(gen);

//     double prob = pow(1 - p, trials);
//     double cum_prob = prob;
//     for (int i = 0; i <= trials; ++i)
//     {
//         if (i >= weight_cap)
//             return weight_cap;
//         if (cum_prob >= u - max_allowed_deviation)
//             return i;
        
//         prob *= (double)(trials - i) / (i + 1) * p / (1 - p);
//         cum_prob += prob;
//     }
// }

//without replacement
pvector<pvector<int> *> sampling(int number_of_trees, pvector<pvector<int> *> &trees, double sum_of_weight, pvector<double> &tree_weights, default_random_engine g)
{
    if (trees.size() < (size_t)number_of_trees)
        return pvector<pvector<int> *>(trees.begin(), trees.end());

    // discrete_distribution<int> dist(tree_weights.begin(), tree_weights.end());
    pvector<pvector<int> *> sampled_trees;
    pvector<bool> visited(trees.size(), 0);
    sampled_trees.reserve(number_of_trees);
    for (int i = 0; i < number_of_trees; i++)
    {
        discrete_distribution<int> dist(tree_weights.begin(), tree_weights.end());
        int tree = dist(g);
        sampled_trees.push_back(trees[tree]);
        visited[tree] = 1;
        tree_weights[tree] = 0;
    }
    for (size_t i = 0; i < visited.size(); i++)
        if (!visited[i])
            delete trees[i];

    return sampled_trees;
}

//if choise = true, estimate with summed weight of a vertex, else minimum weight of a maximum spanning tree
int getUpperbound (const WGraph &g, const pvector<WEdge> &G)
{
    int res = numeric_limits<int>::max();
    int tmp;
    for(auto i : g.vertices())
    {
        tmp = 0;
        for (auto j : g.out_neigh(i))
            tmp += j.w;
        res = min(res,tmp);
    }
    
   
    int n = g.num_nodes();
    pvector<WEdge> t = Kruskal(G, n, false);
    tmp = t[n - 2].v.w;
    tmp *= (n * n);
    res = min(res,tmp);
    
    
    return res;
}

pvector<pvector<WEdge> *> SpanningTreesGenerator(const WGraph &g, const pvector<WEdge> &G, double d, double eps1, double eps2, int n, int m)
{
    assert(eps1 > 0.0 && eps2 > 0.0);
    assert((1.0 - eps2) / (1.0 + eps1) > 2.0 / 3.0);
    double f = 3.0 / 2.0 - (1.0 + eps1) / (1.0 - eps2);
    assert(f > 0);
    double b = ((2 + eps1) * (d + 2.0) * log(n)) / (eps1 * eps1);
    int weight_cap = ceil((1.0 + eps1) * 2.0 * b);
    int c_dash;
    bool lastrun = 0;
    default_random_engine gen;
    // default_random_engine gen(time(NULL));
    const double max_allowed_deviation = 1e-10;
    //Upper bound approximation for mincut value
    c_dash = getUpperbound(g, G);
    
    //console
    cout << c_dash << ", ";
    cout << b * (1 + eps1) << ", \"";
    //console

    while (true)
    {
        //WEdge
        pvector<WEdge> H;
        //remaining_capacity[i] keeps track of edge id i its turn number and remaining capacity
        pvector<pair<int, int>> remaining_capacity;
        //edge_id[i] = corresponding edge id of G in H
        pvector<int> edge_id;
        int M = 0;
        int weight;
        H.reserve(m);
        remaining_capacity.reserve(m);
        edge_id.reserve(m);

        double p = b / c_dash;
        
        for (int i = 0; i < m; i++)
        {
            binomial_distribution<int> dist(G[i].v.w, p);
            weight = dist(gen);
            if(weight > weight_cap)
                weight = weight_cap;
            // int weight = binomial(G[i].v.w, weight_cap, p, gen);
            if (weight != 0)
            {
                M += weight;
                edge_id.push_back(i);
                H.push_back(WEdge(G[i].u, WNode(G[i].v.v, weight)));
                remaining_capacity.push_back(make_pair(weight, weight));
            }
        }

        bool connected;
        {
            auto tmp = WeightedBuilder::Load_CSR_From_Edgelist(H, true);
            connected = isConnected(tmp);
        }

        if (H.size() != 0 && connected)
        {
            pvector<double> tree_weights;
            tree_weights.reserve(m + n * log(n) * log(n));
            //contains packing value and all trees
            double eps2_dash;
            if (p > 1)
                eps2_dash = 1 - ((1-eps2)/(1+eps1));
            else
                eps2_dash = eps2;
            t1.Start();
            pair<double, pvector<pvector<int> *>> res = PackingWeight(H, remaining_capacity, eps2_dash, n, H.size(),M, tree_weights);
            t1.Stop();
            //console
            cout << H.size() << " " << p << " " << t1.Seconds() << " " << res.first << " " << res.second.size() << " ; ";
            //console
            if (res.first < b * (1 + eps1))
                c_dash /= 2.0;
            else
                lastrun = true;

            if (lastrun || p > 1 - max_allowed_deviation)
            {
                int number_of_trees = ceil(-d * log(n) / log(1 - f));
                // pvector<pvector<int> *> tmp;
                // tmp = sampling(number_of_trees, res.second, res.first, tree_weights, gen);
                pvector<pvector<int> *> tmp(res.second.begin(), res.second.end());
                //console
                cout << "\", " << tmp.size() << ", ";
                //console
                pvector<pvector<WEdge> *> trees;
                trees.reserve(tmp.size());
                for (auto i : tmp)
                {
                    pvector<WEdge> *edges = new pvector<WEdge>(n - 1);
                    int k = 0;

                    for (auto idx : *i)
                    {
                        assert(idx < (int)G.size());
                        (*edges)[k] = G[edge_id[idx]];
                        k++;
                    }
                    trees.push_back(edges);
                }
                return trees;
            }
        }
        else
            c_dash /= 2.0;
    }
}

size_t MinCut(const WGraph &g)
{
    t2.Start();
    assert(g.num_edges() > 0);
    int m = g.num_edges();
    int n = g.num_nodes();
    pvector<WEdge> G(m);
    int j = 0;
    for (NodeID u : g.vertices())
    {
        for (WNode wn : g.out_neigh(u))
        {
            assert(wn.w > 0);
            //since g.vertices and g.out_neigh give us sorted NodeID, with u < wn.v we don't count edges twice
            if (u < wn.v)
            {
                G[j] = WEdge(u, wn);
                j++;
            }
        }
    }

    double eps1 = 1.0 / 6.0;
    double eps2 = 1.0 / 5.0;

    //console
    cout << n << "," << m << ",";
    //console

    t.Start();
    pvector<pvector<WEdge> *> tmp = SpanningTreesGenerator(g, G, 1.0, eps1, eps2, n, m);
    t.Stop();

    //console
    cout << t.Seconds() << ", \"";
    // cout << t.Seconds() << " ";
    //console   

    pvector<SpanTree *> trees;
    trees.reserve(tmp.size());
    // vector of csr_trees to keep them in scope
    pvector<WGraph> csr_array(tmp.size());
    j = 0;
    for (auto it : tmp)
    {
        csr_array[j] = WeightedBuilder::Load_CSR_From_Edgelist(*it, true);
        auto xd = new SpanTree(g, csr_array[j++]);
        trees.push_back(xd);
        delete it;
    }

    t.Start();
    int res = numeric_limits<int>::max();
    // 1,2,3 indiciates incomparablecut, comparablecut or 1-respect mincut
    int which_solution = 0;
    for (auto i : trees)
    {
        int tmp = i->compute(which_solution);
        //console
        // cout << tmp << " ";
        //console
        if (res > tmp)
            res = tmp;
    }
    
    t.Stop();
    // PrintStep("min compute:                                               ", t.Seconds());
    for (auto i : trees)
        delete i;
    // cout << endl <<"res: " << res << " solution of: " << which_solution << endl;
    t2.Stop();

    //console
    cout << "\", "<< t.Seconds() << ", " << t2.Seconds() << ", " << "\"" << res << " ";
    // console

    return res;
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
    // if (mincut != test_min)
    //     cout << mincut << " != " << test_min << endl;
    cout << mincut << "\"\n";
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
    if (!isConnected(g))
    {
        cout << "Input graph is not connected" << endl;
        return -2;
    }

    //console
    cout << "nodes, edges, mincut estimate, packing treshold, \"H size, p, packing time, packing value, packing size\", sample size, trees generator time, cut values of sampled trees, mincut computation time, overall time, end result vs correct result\n";
    // cout << "nodes, edges, mincut estimate, packing treshold, \"H size, p, packing time, packing value, packing size\", sample size, trees generator time\n";
    //console

    BenchmarkKernel(cli, g, MinCut, PrintMinCutValue, MINCUTVerifier);
}
