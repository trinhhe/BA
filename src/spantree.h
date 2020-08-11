#ifndef SPANTREE_H_
#define SPANTREE_H_

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <string.h>
#include <assert.h>
#include <queue>
#include <stack>
#include <unordered_map>
#include <set>
#include <list>
// #include <typeinfo>

#include "graph.h"
#include "pvector.h"
#include "minimumpath.h"
#include "timer.h"

/*
Author: Henry Trinh

class to represent spanntree
*/
using namespace std;

typedef int32_t NodeID;
typedef int32_t WeightT;
typedef NodeWeight<NodeID, WeightT> WNode;
typedef EdgePair<NodeID, WNode> WEdge;
typedef CSRGraph<NodeID, WNode> WGraph;
typedef int32_t PathID;
int INT_MAX = numeric_limits<int>::max();
int INT_MIN = numeric_limits<int>::min();

class SpanTree
{
    const WGraph &graph;
    const WGraph &tree;
    int root;
    int num_nodes;
    int num_paths;
    //stores parent of a node, parent of root is -1
    pvector<NodeID> parent;
    //contains path ID where NodeID belongs to (maybe not needed)
    pvector<NodeID> path_id;
    //path[i] = j means that j is next node of i in some path, if i is last node in chain path[i] = -1
    pvector<NodeID> path;
    //head[i] is the head vertex of a path (head vertex is closer to root of tree)
    pvector<NodeID> head;
    //stores position of a node in MinPath structure
    pvector<NodeID> minpathPos;
    //length of paths
    pvector<PathID> lengths;

    //pointers to all minpath structures of the spantree
    pvector<MinimumPath *> P;

    //vector of 1-min cut values of all nodes
    vector<int> C;
    //rho[i] holds weight of all edges of the descandants where both endpoints share i as lca
    vector<int> rho;
    //help vectors for LCA queries
    pvector<NodeID> levels;
    pvector<NodeID> inlabel;
    pvector<NodeID> ascendant;
    unordered_map<int, int> head_;
    //postorder[i] stores the node visited on postorder iteration i (needed for computing treefix sum)
    pvector<NodeID> postorder;

public:
    // CONSTRUCTORS //

    SpanTree(const WGraph &G, const WGraph &T, int root, int num_nodes, int num_paths, pvector<NodeID> &parent, 
        pvector<NodeID> &path_id, pvector<NodeID> &path, pvector<NodeID> &head, pvector<NodeID> &minpathPos, pvector<PathID> &lengths) : 
        graph(G), tree(T), root(root),num_nodes(num_nodes), num_paths(num_paths), parent(parent.begin(), parent.end()), path_id(path_id.begin(), path_id.end()), path(path.begin(), 
        path.end()), head(head.begin(), head.end()), minpathPos(minpathPos.begin(), minpathPos.end()), lengths(lengths.begin(), lengths.end()), C(num_nodes), rho(num_nodes, 0), levels(num_nodes), inlabel(num_nodes), ascendant(num_nodes), postorder(num_nodes)
    {
        P.reserve(num_paths);
        for (int i = 0; i < num_paths; i++)
        {
            P.push_back(new MinimumPath(lengths[i]));
        }
    }

    SpanTree(const WGraph &G, const WGraph &T) : graph(G), tree(T), num_nodes(G.num_nodes()), num_paths(0), parent(num_nodes, -1),
                                                 path_id(num_nodes, -1), path(num_nodes), head(num_nodes), minpathPos(num_nodes),
                                                 C(num_nodes), rho(num_nodes, 0), levels(num_nodes), inlabel(num_nodes), ascendant(num_nodes), postorder(num_nodes)
    {
        PathSegmentation();
        P.reserve(num_paths);
        for (int i = 0; i < num_paths; i++)
        {
            P.push_back(new MinimumPath(lengths[i]));
        }
    }


    ~SpanTree()
    {
        for (auto p : P)
            delete p;
        P.clear();
        // cout << "Destructor called" << endl;
    }

    void print()
    {
        cout << "\n------------------------------------------" << endl;
        cout << "Graph Topology: " << endl;
        graph.PrintTopology();
        cout << "Tree Topology: " << endl;
        tree.PrintTopology();
        cout << "        parent, path, path_id, head, minpathPos: " << endl;
        for (int i = 0; i < num_nodes; i++)
        {
            cout << i << "      " << parent[i] << "      " << path[i] << "      " << path_id[i] << "       " << head[i] << "       " << minpathPos[i] << endl;
        }
        cout << "num_paths: " << num_paths << endl;
        cout << "root: " << root << endl;
        // cout << "P: " << endl;
        // int j = 0;
        // for (auto i : P)
        // {
        //     cout << "Path ID: " << j++ << endl;
        //     i->print();
        //     cout << endl;
        // }
        cout << "------------------------------------------" << endl;
    }
    void printpath()
    {
        cout << "\n------------------------------------------" << endl;
        int j = 0;
        for (auto i : P)
        {
            cout << "Path ID: " << j++ << endl;
            i->print();
            cout << endl;
        }
        cout << "------------------------------------------" << endl;
    }

    //function to trace a path until we hit a branching node
    void path_delving(queue<NodeID> &q, pvector<NodeID> &vertex_degree, pvector<NodeID> &temp_vertex_degree,
                      pvector<bool> &visited, NodeID currentNode, NodeID u)
    {

        path_id[u] = path_id[currentNode];
        path[u] = currentNode;
        visited[u] = true;

        for (NodeID v : tree.out_neigh(u))
        {
            if (!visited[v])
            {
                parent[u] = v;
                //parent is non-branching node
                if (vertex_degree[v] == 2)
                {
                    path_delving(q, vertex_degree, temp_vertex_degree, visited, u, v);
                }
                //parent has multiple children
                else if (vertex_degree[v] > 2)
                {
                    temp_vertex_degree[v]--;
                    if (temp_vertex_degree[v] == 1)
                        //change to emplace instead of push
                        q.push(v);
                }
                else
                {
                    path_id[v] = path_id[u];
                    path[v] = u;
                    visited[v] = true;
                }
            }
        }
    }
    //start from leaves and go up until we hit root
    void PathSegmentation()
    {
        //storage of degrees for every bought phase
        pvector<NodeID> vertex_degree(num_nodes);
        //temporary storage of degrees which keeps track after! bought phases
        pvector<NodeID> temp_vertex_degree(num_nodes);
        pvector<bool> visited(num_nodes, false);
        queue<NodeID> q;
        for (NodeID i : tree.vertices())
        {
            vertex_degree[i] = tree.out_degree(i);
            temp_vertex_degree[i] = tree.out_degree(i);
            //fill queue with all leaves
            if (tree.out_degree(i) == 1)
            {
                q.push(i);
            }
        }
        //path ID counter
        PathID j = 0;
        int queue_size = q.size();
        while (queue_size > 0)
        {
            for (int i = 0; i < queue_size; ++i)
            {
                NodeID currentLeaf = q.front();
                q.pop();
                if (!visited[currentLeaf])
                {
                    path_id[currentLeaf] = j++;
                    path[currentLeaf] = -1;
                    visited[currentLeaf] = true;
                }

                for (NodeID u : tree.out_neigh(currentLeaf))
                {
                    if (!visited[u])
                    {
                        parent[currentLeaf] = u;
                        //parent is non-branching node
                        if (vertex_degree[u] == 2)
                        {
                            path_delving(q, vertex_degree, temp_vertex_degree, visited, currentLeaf, u);
                        }
                        //parent has multiple children
                        else if (vertex_degree[u] > 2)
                        {
                            temp_vertex_degree[u]--;
                            if (temp_vertex_degree[u] == 1)
                                q.push(u);
                        }
                        //only a path left
                        else
                        {
                            path_id[u] = path_id[currentLeaf];
                            path[u] = currentLeaf;
                            visited[u] = true;
                        }
                    }
                }
            }
            //we've gone through all leaves in this bough phase, enter to next bought phase

            queue_size = q.size();
            for (int k = 0; k < num_nodes; ++k)
            {
                vertex_degree[k] = temp_vertex_degree[k];
            }
        }

        // determine head vertex, root, num_paths and minpathPos of every vertex
        for (NodeID i = 0; i < num_nodes; i++)
        {
            if (parent[i] == -1)
                root = i;
            if (path_id[i] > num_paths)
                num_paths = path_id[i];
            //if i is root or start vertex of a path
            if (parent[i] == -1 || path[parent[i]] != i)
            {
                for (NodeID j = i, currentPos = 0; j != -1; j = path[j])
                {
                    head[j] = i;
                    minpathPos[j] = currentPos++;
                }
            }
        }
        num_paths++;
        lengths.reserve(num_paths);
        for (NodeID i = 0; i < num_nodes; i++)
        {
            if (path[i] == -1)
            {
                lengths[path_id[i]] = minpathPos[i] + 1;
            }
        }
    }

    //build on existing leave weights the Minpath structures
    void build()
    {
        for (auto &i : P)
        {
            i->build();
        }
    }
    // function that will fill the minpath leaves with the corresponding initial weight (Karger, smallest cut that 1-respects given tree)
    void set(NodeID v, const int &value)
    {
        P[path_id[v]]->set(minpathPos[v], value);
    }

    /*weight initialization 
    See "Minimum Cuts in Near Linear Time" by David R. Karger 
    + "On Finding Lowest Common Ancestors: Simplification and Parallelization" by Baruch Schieber and Uzi Vishkin
    */
    void InitializeWeight()
    {
        //preprocessing O(n) to query LCA in O(1)
        pvector<bool> visited(num_nodes, false);
        pvector<int> preorder(num_nodes);
        pvector<int> sizes(num_nodes);
        
        int index = 1;
        PreOrder(root, index, 0, preorder, levels, sizes, visited);

        //compute inlabel
        for (NodeID v = 0; v < num_nodes; v++)
        {
            //step  2.1
            int i = preorder[v] - 1;
            i = i ^ (preorder[v] + sizes[v] - 1);
            i = log2(i);
            //step 2.2.
            inlabel[v] = ((preorder[v] + sizes[v] - 1) >> i) << i;
        }
        //compute ascendant
        GetAscendants(root, ascendant, inlabel);

        //compute head
        for (NodeID v = 0; v < num_nodes; v++)
        {
            if (inlabel[v] != inlabel[parent[v]] || parent[v] == -1)
                head_[inlabel[v]] = v;
        }
        //preprocessing done for lca

        //delta[i] holds weight of all edges of the descandants
        vector<NodeID> delta(num_nodes, 0);
        PostOrder(root, postorder);

        for (NodeID v : graph.vertices())
            for (WNode u : graph.out_neigh(v))
                delta[v] += u.w;

        for (NodeID v : graph.vertices())
            for (WNode u : graph.out_neigh(v))
            {
                if (v < u.v)
                {
                    NodeID z = LCA_Query(u.v, v, inlabel, ascendant, levels, head_);
                    // cout << v << " " << u.v << " " << z << endl;
                    rho[z] += u.w;
                }
            }

        for (NodeID v : postorder)
        {
            if (v != root)
            {
                delta[parent[v]] += delta[v];
                rho[parent[v]] += rho[v];
            }
        }

        for (int i = 0; i < num_nodes; i++)
        {
            set(i, delta[i] - 2 * rho[i]);
            C[i] = delta[i] - 2 * rho[i];
        }
        build();
        // int j = 0;
        // cout << "preorder: " << endl;
        // for (auto i : preorder)
        //     cout << j++ << ": " << i << endl;

        // j = 0;
        // cout << "postorder: " << endl;
        // for (auto i : postorder)
        //     cout << j++ << ": " << i << endl;

        // j = 0;
        // cout << "delta: " << endl;
        // for (auto i : delta)
        //     cout << j++ << ": " << i << endl;
        // j = 0;
        // cout << "rho: " << endl;
        // for (auto i : rho)
        //     cout << j++ << ": " << i << endl;
        // j = 0;
        // cout << "weights: " << endl;
        // for (int i = 0; i < num_nodes; i++)
        //     cout << j++ << ": " << C[i] << endl;
    }

    /* ------------------------------functions for lca----------------------------*/

    //build preorder numbering, size of subtrees and depth level of all vertices
    void PreOrder(int u, int &index, int level, pvector<int> &preorder, pvector<int> &levels, pvector<int> &sizes, pvector<bool> &visited)
    {
        if (!visited[u])
            visited[u] = true;
        levels[u] = level;
        preorder[u] = index++;
        sizes[u] = 1;

        for (auto it : tree.out_neigh(u))
        {
            if (!visited[it])
            {
                PreOrder(it, index, level + 1, preorder, levels, sizes, visited);
                sizes[u] += sizes[it];
            }
        }

        //iterative but can't get subtree sizes

        // pvector<int> visited(num_nodes, -1);
        // int index = 1;
        // stack<int> vertices;
        // vertices.push(root);
        // while (!vertices.empty())
        // {
        //     int curr = vertices.top();
        //     vertices.pop();
        //     if (visited[curr] == -1)
        //     {
        //         visited[curr] = 0;
        //         preorder[curr] = index++;
        //         for (auto it : tree.out_neigh(curr))
        //             vertices.push(it);
        //     }
        // }
    }

    //compute ascendants like in paper in bfs traversal
    void GetAscendants(const int &root, pvector<int> &ascendant, pvector<int> &inlabel)
    {
        queue<int> vertices;
        pvector<bool> visited(num_nodes, false);
        visited[root] = true;
        vertices.push(root);
        while (!vertices.empty())
        {
            int curr = vertices.front();
            vertices.pop();
            //compute ascendant
            if (curr == root)
            {
                int l = ceil(log2(num_nodes + 1)) - 1;
                ascendant[root] = 1 << l;
            }
            else if (inlabel[curr] == inlabel[parent[curr]])
                ascendant[curr] = ascendant[parent[curr]];
            else
            {
                int i = log2(inlabel[curr] - (inlabel[curr] & (inlabel[curr] - 1)));
                ascendant[curr] = ascendant[parent[curr]] + (1 << i);
            }
            for (auto it : tree.out_neigh(curr))
            {
                if (!visited[it])
                {
                    visited[it] = true;
                    vertices.push(it);
                }
            }
        }
    }

    //O(1) lca query
    NodeID LCA_Query(const int &u, const int &v, const pvector<int> &inlabel,
                     const pvector<int> &ascendant, const pvector<int> &level, const unordered_map<int, int> &head)
    {
        assert(u >= 0 && v >= 0 && u < num_nodes && v < num_nodes);
        //Case A
        if (inlabel[u] == inlabel[v])
            return level[u] <= level[v] ? u : v;
        //Case B
        else
        {
            int i, common, j, inlabel_z, k, u_hat, v_hat, inlabel_w, w;
            //step 1
            i = log2(inlabel[u] ^ inlabel[v]);
            //dunno why this is in the paper in step 1
            // b = ((inlabel[u] >> (i + 1)) << (i + 1)) + (1 << i);

            //step 2
            common = ascendant[u] & ascendant[v];
            common = (common >> i) << i;
            j = log2(common - (common & (common - 1)));
            inlabel_z = ((inlabel[u] >> (j + 1)) << (j + 1)) + (1 << j);
            //step 3
            if (inlabel[u] == inlabel_z)
                u_hat = u;
            else
            {
                k = ((1 << j) - 1) & ascendant[u];
                if (k != 0)
                    k = log2(k);
                inlabel_w = ((inlabel[u] >> (k + 1)) << (k + 1)) + (1 << k);
                w = head.find(inlabel_w)->second;
                u_hat = parent[w];
            }

            if (inlabel[v] == inlabel_z)
                v_hat = v;
            else
            {
                k = ((1 << j) - 1) & ascendant[v];
                if (k != 0)
                    k = log2(k);
                inlabel_w = ((inlabel[v] >> (k + 1)) << (k + 1)) + (1 << k);
                w = head.find(inlabel_w)->second;
                v_hat = parent[w];
            }

            return level[u_hat] <= level[v_hat] ? u_hat : v_hat;
        }
    }

    /* --------------------------end of functions for lca----------------------------*/

    /* ----------------------functions for treefix sum (needed for 1-mincut values)----------------------------*/

    void PostOrder(const NodeID &root, pvector<NodeID> &postorder)
    {
        pvector<bool> visited(num_nodes, false);
        int index = 0;
        stack<NodeID> st;
        st.push(root);
        while (!st.empty())
        {
            NodeID current = st.top();
            if (!visited[current])
            {
                visited[current] = true;
                for (NodeID i : tree.out_neigh(current))
                    if (!visited[i])
                        st.push(i);
            }
            else
            {
                postorder[index++] = current;
                st.pop();
            }
        }
    }

    /* -----------------------end of functions for treefix sum (needed for 1-mincut values)-----------------------*/

    int MinPath(NodeID v)
    {
        int res = INT_MAX;
        //go through all paths until we hit the root of spantree
        for (; head[v] != root; v = parent[head[v]])
        {
            res = min(res, P[path_id[v]]->minprefix(minpathPos[v]));
        }
        res = min(res, P[path_id[v]]->minprefix(minpathPos[v]));
        return res;
    }

    void AddPath(NodeID v, const int &value)
    {
        //go through all paths until we hit the root of spantree
        for (; head[v] != root; v = parent[head[v]])
        {
            P[path_id[v]]->addprefix(minpathPos[v], value);
        }
        P[path_id[v]]->addprefix(minpathPos[v], value);
    }

    // 1.Case (u,v) and (s,t) cuts are incomparable. (u is not ancestor of s or vice versa)
    int IncomparableCut()
    {
        //adj list for now, csr format with its continious memory layout can't handle a case where a vertex has more edges after a contraction
        list<WNode> adj[num_nodes];
        for (NodeID u : graph.vertices())
        {
            for (WNode wn : graph.out_neigh(u))
            {
                adj[u].push_back(wn);
            }
        }
        list<WNode>::iterator it1;
        list<WNode>::iterator it2;

        int pos_inf = (int)2.5e8;
        //store smallest value from MinPath queries
        pvector<NodeID> minvalues(num_nodes, INT_MAX);
        int tmp_min;
        //visited[i] = NodeID i is true indicates to to ignore all edges which have i as a node (contract edges)
        pvector<NodeID> visited(num_nodes, 0);
        //keep track which nodes to contract after bough-phase
        queue<NodeID> to_remove;
        //keep track of vertex degrees to decide the leaves in next bough phase
        vector<NodeID> degree(num_nodes);
        // leaves in current bough phase
        std::set<NodeID> leaves;
        for (NodeID i = 0; i < num_nodes; i++)
        {
            degree[i] = tree.out_degree(i);
            //in case tree is a path, don't include root as a leaf
            if (tree.out_degree(i) == 1 && i != root)
                leaves.insert(i);
        }

        int iterator = num_paths;
        while (iterator > 0)
        {
            for (NodeID v : leaves)
            {
                //keep track of the order of vertex in a bough to revert the AddPath operations
                stack<pair<NodeID, int>> stack;
                stack.push(make_pair(v, pos_inf));
                AddPath(v, pos_inf);

                for (; head[v] != v; v = parent[v])
                {
                    degree[v] -= 2;
                    visited[v] = true;
                    to_remove.push(v);
                    for (WNode x : adj[v])
                    {
                        stack.push(make_pair(x, -2 * x.w));
                        AddPath(x, -2 * x.w);
                    }
                    for (WNode x : adj[v])
                    {
                        tmp_min = MinPath(x);
                        if (minvalues[v] > tmp_min)
                            minvalues[v] = tmp_min;
                    }
                }

                //include head vertex of bough too
                degree[v] -= 2;
                visited[v] = true;
                to_remove.push(v);
                for (WNode x : adj[v])
                {
                    stack.push(make_pair(x, -2 * x.w));
                    AddPath(x, -2 * x.w);
                }
                for (WNode x : adj[v])
                {
                    tmp_min = MinPath(x);
                    if (minvalues[v] > tmp_min)
                        minvalues[v] = tmp_min;
                }
                // "remove" edge to head of vertex, later on all nodes with degree = 1 will be new leaves
                if (v != root)
                    degree[parent[v]]--;
                //reverse operations
                while (!stack.empty())
                {
                    pair<NodeID, int> curr = stack.top();
                    stack.pop();
                    AddPath(curr.first, -1 * curr.second);
                }

                iterator--;
            }
            leaves.clear();

            while (!to_remove.empty())
            {
                NodeID curr = to_remove.front();
                to_remove.pop();
                NodeID x = parent[head[curr]];
                for (it1 = adj[curr].begin(); it1 != adj[curr].end(); ++it1)
                {
                    NodeID y = it1->v;
                    //if visited[y] = true both endpoints of an edge were in a bough,
                    if (visited[y])
                    {
                        //if endpoints are not in the same bough we "contract" the edge to the parent of the bough.
                        //we don't add the other edge to have "symmetry" in adj list since the other endpoint will also be processed later on
                        y = parent[head[y]];
                        if (x != y)
                            adj[x].push_back(WNode(y, it1->w));
                    }
                    else
                    {
                        if (x != y)
                        {
                            adj[x].push_back(WNode(y, it1->w));
                            adj[y].push_back(WNode(x, it1->w));
                        }
                        for (it2 = adj[y].begin(); it2 != adj[y].end(); ++it2)
                            if (it2->v == curr)
                            {
                                adj[y].erase(it2);
                                break;
                            }
                    }
                }
                adj[curr].clear();
                if (head[curr] != root && degree[parent[head[curr]]] == 1 && path[parent[head[curr]]] == -1)
                    leaves.insert(parent[head[curr]]);
            }
            //one node left after contractions which is the root
            if (iterator == 1 && leaves.empty())
                leaves.insert(root);
        }

        int res = INT_MAX;
        for (size_t i = 0; i < minvalues.size(); i++)
        {
            // assert(
            //     ((minvalues[i] > 0) && (C[i] > INT_MAX - minvalues[i])) ||
            //     ((minvalues[i] < 0) && (C[i] < INT_MIN - minvalues[i])));
            if (res > minvalues[i] + C[i])
                res = minvalues[i] + C[i];
        }

        return res;
    }

    // 2.Case: Difference of two cuts
    int ComparableCut()
    {
        list<WNode> adj[num_nodes];
        for (NodeID u : graph.vertices())
        {
            for (WNode wn : graph.out_neigh(u))
            {
                adj[u].push_back(wn);
            }
        }
        list<WNode>::iterator it1;
        list<WNode>::iterator it2;

        int pos_inf = (int)2.5e8;

        //store smallest value from MinPath queries
        vector<NodeID> minvalues(num_nodes, INT_MAX);
        //visited[i] = NodeID i is true indicates to to ignore all edges which have i as a node (contract edges)
        vector<NodeID> visited(num_nodes, 0);
        //keep track which nodes to contract after bough-phase
        queue<NodeID> to_remove;
        //keep track of vertex degrees to decide the leaves in next bough phase
        vector<NodeID> degree(num_nodes);
        // leaves in current bough phase
        std::set<NodeID> leaves;
        for (NodeID i = 0; i < num_nodes; i++)
        {
            degree[i] = tree.out_degree(i);
            //in case tree is a path, don't include root as a leaf
            if (tree.out_degree(i) == 1 && i != root)
                leaves.insert(i);
        }
        int iterator = num_paths;
        //do not consider C(root) 
        AddPath(root, pos_inf);
        while (iterator > 0)
        {
            for (NodeID v : leaves)
            {
                //keep track of the order of vertex in a bough to revert the AddPath operations
                stack<pair<NodeID, int>> stack;

                for (; head[v] != v; v = parent[v])
                {
                    degree[v] -= 2;
                    visited[v] = true;
                    to_remove.push(v);
                    for (WNode x : adj[v])
                    {
                        stack.push(make_pair(x, 2 * x.w));
                        AddPath(x, 2 * x.w);
                    }
                    minvalues[v] = MinPath(parent[v]);
                }

                //include head vertex of bough too
                degree[v] -= 2;
                visited[v] = true;
                to_remove.push(v);
                for (WNode x : adj[v])
                {
                    stack.push(make_pair(x, 2 * x.w));
                    AddPath(x, 2 * x.w);
                }
                if(v != root){
                    minvalues[v] = MinPath(parent[v]);
                    // "remove" edge to head of vertex, later on all nodes with degree = 1 will be new leaves
                    degree[parent[v]]--;
                }
                    
                //reverse operations
                while (!stack.empty())
                {
                    pair<NodeID, int> curr = stack.top();
                    stack.pop();
                    AddPath(curr.first, -1 * curr.second);
                }

                iterator--;
            }
            leaves.clear();

            while (!to_remove.empty())
            {
                NodeID curr = to_remove.front();
                to_remove.pop();
                NodeID x = parent[head[curr]];
                for (it1 = adj[curr].begin(); it1 != adj[curr].end(); ++it1)
                {
                    NodeID y = it1->v;
                    //if visited[y] = true both endpoints of an edge were in a bough,
                    if (visited[y])
                    {
                        //if endpoints are not in the same bough we "contract" the edge to the parent of the bough.
                        //we don't add the other edge to have "symmetry" in adj list since the other endpoint will also be processed later on
                        y = parent[head[y]];
                        if (x != y)
                            adj[x].push_back(WNode(y, it1->w));
                    }
                    else
                    {
                        if (x != y)
                        {
                            adj[x].push_back(WNode(y, it1->w));
                            adj[y].push_back(WNode(x, it1->w));
                        }
                        for (it2 = adj[y].begin(); it2 != adj[y].end(); ++it2)
                            if (it2->v == curr)
                            {
                                adj[y].erase(it2);
                                break;
                            }
                    }
                }
                adj[curr].clear();
                if (head[curr] != root && degree[parent[head[curr]]] == 1 && path[parent[head[curr]]] == -1)
                    leaves.insert(parent[head[curr]]);
            }
            //one node left after contractions which is the root
            if (iterator == 1 && leaves.empty())
                leaves.insert(root);
            
            //recompute rho values for nodes after contraction which are still needed
            for (int i = 0; i < num_nodes; i++)
                if(!visited[i])
                    rho[i] = 0;
            for (int v = 0; v < num_nodes; v++)
                for (it1 = adj[v].begin(); it1 != adj[v].end(); it1++)
                    if (v < it1->v)
                    {
                        NodeID z = LCA_Query(it1->v, v, inlabel, ascendant, levels, head_);
                        rho[z] += it1->w;
                    }
            for (NodeID v : postorder)
                if (v != root && !visited[v])
                    rho[parent[v]] += rho[v];
        }

        int res = INT_MAX;
        // cout << "\nminvalues: " << endl;
        // for (auto i : minvalues)
        //     cout << i << " ";
        // cout << "\nrho: " << endl;
        // for (auto i : rho)
        //     cout << i << " ";
        // cout << "\nC: " << endl;
        // for (auto i : C)
        //     cout << i << " ";

        for (size_t i = 0; i < minvalues.size(); i++)
        {
            // assert(
            //     ((minvalues[i] > 0) && (C[i] > INT_MAX - minvalues[i])) ||
            //     ((minvalues[i] < 0) && (C[i] < INT_MIN - minvalues[i])));
            if (res > minvalues[i] - 4 * rho[i] - C[i])
                res = minvalues[i] - 4 * rho[i] - C[i];
        }

        return res;
    }

    int compute()
    {
        Timer t;
        if (num_nodes == 1)
            return 0;

        int res = INT_MAX;
        t.Start();
        InitializeWeight();
        t.Stop();
        // print();
        // PrintStep('InitializeWeight', t.Seconds());
        int c1 = IncomparableCut();
        int c2 = ComparableCut();
        int tmp;
        if (c1 > 0 && c2 > 0)
            tmp = min(c1, c2);
        else if(c1 > 0 && c2 < 0)
            tmp = c1;
        else
            tmp = c2;
        if (res > tmp && tmp > 0)
            res = tmp;
        // int lol = INT_MAX;
        // int index;
        // cout << "C(V): ";
        for (int i = 0; i < num_nodes; i++)
        {
            // cout << C[i] << " ";
            if (C[i] < res && C[i] != 0)
                res = C[i];
            // if (C[i] < lol && C[i] != 0)
            //     lol = C[i], index = i;
        }
        //console
        // cout << "values " <<lol << " " << index << " ";
        // cout << c1 << " " << c2 << " " << res << endl;
        //console

        // print();
        return res;
    }
};

#endif