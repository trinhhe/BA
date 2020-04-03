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
    WGraph &g;
    int num_nodes;
    int num_paths;
    //stores parent of a node, parent of root is -1
    pvector<NodeID> parent;
    //contains path ID where NodeID belongs to
    pvector<NodeID> path_id;
    // pvector<PathID> parent_path;
    // pvector<PathID> child_path;

public:
    // CONSTRUCTORS //

    SpanTree(WGraph &g) : g(g), num_nodes(g.num_nodes()), parent(num_nodes, -1), path_id(num_nodes, -1)
    {
        PathSegmentation();
    }

    ~SpanTree() { cout << "Destructor called" << endl; }

    void print()
    {
        cout << "Tree Topology: " << endl;
        g.PrintTopology();
        cout << "parent array: " << endl;
        for (int i = 0; i < num_nodes; i++)
        {
            cout << i << " " << parent[i] << " " << path_id[i] << endl;
        }
    }

    void PathSegmentation()
    {
        //temporary storage of degrees
        pvector<NodeID> vertex_degree(num_nodes);
        //helper vector to follow down nodes on same path id
        pvector<NodeID> childpath(num_nodes, -1);
        queue<NodeID> q;
        int v = num_nodes;
        for (NodeID i : g.vertices())
        {
            vertex_degree[i] = g.out_degree(i);
            //fill queue with all leaves
            if (g.out_degree(i) == 1)
            {
                q.push(i);
                --v;
            }
        }
        PathID j = -1;
        //loop until total vertex is less than 2
        while (v > 0)
        {
            int queue_size = q.size();
            //process all leaves
            for (int i = 0; i < queue_size; i++)
            {
                NodeID currentLeaf = q.front();
                q.pop();
                // -1 means node hasn't path id yet
                if (path_id[currentLeaf] == -1)
                {
                    j++;
                    path_id[currentLeaf] = j;
                }
                //decrease degree for each neighbour of currentLeaf
                for (NodeID u : g.out_neigh(currentLeaf))
                {
                    //since we go from leaves to root and a node can only have one parent,
                    //the node with degree higher 1 must be the parent
                    if (vertex_degree[u] > 1)
                    {
                        parent[currentLeaf] = u;
                    }

                    if (vertex_degree[u] > 2 && path_id[u] == -1)
                    {
                        j++;
                        path_id[u] = j;
                    }
                    else if (vertex_degree[u] == 2 && path_id[u] == -1)
                    {
                        path_id[u] = j;
                        childpath[u] = currentLeaf;
                        // path_id[u] = path_id[currentLeaf];
                    }
                    else if (vertex_degree[u] == 2 && path_id[u] != -1)
                    {
                        //higher node alrdy has path id
                        //go down and rename all children of the path id
                        childpath[u] = currentLeaf;
                        NodeID k = childpath[u];
                        while (k != -1)
                        {
                            path_id[k] = path_id[u];
                            k = childpath[k];
                        }
                    }

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
    }

    void MinPath()
    {
    }

    void AddPath()
    {
    }
};
#endif