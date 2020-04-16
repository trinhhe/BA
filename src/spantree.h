#ifndef SPANTREE_H_
#define SPANTREE_H_

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <string.h>
#include <assert.h>
#include <queue>

#include "graph.h"
#include "pvector.h"
#include "minimumpath.h"

/*
Author: Henry Trinh

class to represent spanntree
*/
using namespace std;

typedef int32_t NodeID;
typedef int32_t WeightT;
typedef NodeWeight<NodeID, WeightT> WNode;
typedef CSRGraph<NodeID, WNode> WGraph;
typedef int32_t PathID;

class SpanTree
{
    const WGraph &tree;
    int root;
    int num_nodes;
    // int num_paths;
    //stores parent of a node, parent of root is -1
    pvector<NodeID> parent;
    //contains path ID where NodeID belongs to (maybe not needed)
    pvector<NodeID> path_id;
    //path[i] = j means that j is next node of i in some path, if i is last node in chain path[i] = -1
    pvector<NodeID> path;
    //head[i] is the head vertex of a path which i belongs to (head vertex is closer to root of tree)
    pvector<NodeID> head;
    //stores position of a node in MinPath structure
    pvector<NodeID> minpathPos;

    //tree structure for all paths
    MinimumPath P;

public:
    // CONSTRUCTORS //

    SpanTree(WGraph &g) : tree(g), num_nodes(g.num_nodes()), parent(num_nodes, -1),
                          path_id(num_nodes, -1), path(num_nodes), minpathPos(num_nodes), P(num_nodes)
    {
        PathSegmentation();
        // num_paths = head.size();
    }

    ~SpanTree() { cout << "Destructor called" << endl; }

    void print()
    {
        cout << "Tree Topology: " << endl;
        tree.PrintTopology();
        cout << "parent array, path, head and minpathPos : " << endl;
        for (int i = 0; i < num_nodes; i++)
        {
            cout << i << " " << parent[i] << " " << path[i] << " " << head[i] << " " << minpathPos[i] << endl;
        }
    }

    //start from leaves and go up until we hit root
    void PathSegmentation()
    {
        //temporary storage of degrees
        pvector<NodeID> vertex_degree(num_nodes);
        queue<NodeID> q;
        int v = num_nodes;
        PathID j = 0;
        for (NodeID i : tree.vertices())
        {
            vertex_degree[i] = tree.out_degree(i);
            //fill queue with all leaves
            if (tree.out_degree(i) == 1)
            {
                q.push(i);
                path_id[i] = j;
                // tail.push_back(i);
                path[i] = -1;
                head.push_back(i);
                ++j;
                --v;
            }
        }
        NodeID currentLeaf = -1;
        //loop until total vertex is less than 2
        while (v > 0)
        {
            int queue_size = q.size();
            //process all leaves
            for (int i = 0; i < queue_size; i++)
            {
                currentLeaf = q.front();
                q.pop();

                for (NodeID u : tree.out_neigh(currentLeaf))
                {
                    //since we go from leaves to root and a node can only have one parent,
                    //the node with degree higher 1 must be the parent
                    if (vertex_degree[u] > 1)
                    {
                        parent[currentLeaf] = u;
                    }

                    if (vertex_degree[u] == 2 && path_id[u] == -1)
                    {
                        path_id[u] = path_id[currentLeaf];
                        head[path_id[currentLeaf]] = u;
                        path[u] = currentLeaf;
                    }
                    //decrease degree for each neighbour of currentLeaf
                    vertex_degree[u]--;
                    //insert into queue if node becomes leaf
                    if (vertex_degree[u] == 1)
                    {
                        q.push(u);
                        --v;
                    }
                }
            }
        }

        //two possible candidates for tree root
        if (q.size() == 2)
        {
            int child = q.front();
            q.pop();
            parent[child] = q.front();
        }
        root = q.front();
        q.pop();

        // add parent to last node in iteration from before
        for (NodeID u : tree.out_neigh(root))
            if (parent[u] == -1)
                parent[u] = root;

        //determine head vertex and minpathPos of every vertex
        for (NodeID i = 0, currentPos = 0; i < num_nodes; i++)
        {
            //if i is root or start vertex of a path
            if (parent[i] == -1 || path[parent[i]] != i)
                for (NodeID j = i; j != -1; j = path[j])
                {
                    head[j] = i;
                    minpathPos[j] = currentPos++;
                }
        }
    }
    //function that will fill the minpath leaves with the corresponding initial weight (Karger, smallest cut that 1-respects given tree)
    void set(NodeID v, const int &value)
    {
        P.set(minpathPos[v], value);
    }

    int MinPath(NodeID v)
    {
        int res = numeric_limits<int>::max();
        //go through all paths until we hit the root of spantree
        for (; head[v] != root; v = parent[head[v]])
        {
            res = min(res, P.minprefix(minpathPos[v]));
        }
        res = min(res, P.minprefix(minpathPos[v]));
        return res;
    }

    void AddPath(NodeID v, const int &value)
    {
        //go through all paths until we hit the root of spantree
        for (; head[v] != root; v = parent[head[v]])
        {
            P.addprefix(minpathPos[v], value);
        }
        P.addprefix(minpathPos[v], value);
    }
};
#endif