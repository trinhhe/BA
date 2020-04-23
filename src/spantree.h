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

    //tree structure for all paths
    pvector<MinimumPath> P;

public:
    // CONSTRUCTORS //

    SpanTree(WGraph &g) : tree(g), num_nodes(g.num_nodes()), num_paths(0), parent(num_nodes, -1),
                          path_id(num_nodes, -1), path(num_nodes), head(num_nodes), minpathPos(num_nodes)
    {
        PathSegmentation();
        // num_paths = head.size();
    }

    ~SpanTree() { cout << "Destructor called" << endl; }

    void print()
    {
        cout << "Tree Topology: " << endl;
        tree.PrintTopology();
        cout << "parent, path, path_id, head: " << endl;
        for (int i = 0; i < num_nodes; i++)
        {
            cout << i << " " << parent[i] << " " << path[i] << " " << path_id[i] << " " << head[i] << endl;
        }
        cout << "num_paths: " << num_paths << endl;
        cout << "root: " << root << endl;
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
            //if i is root or start vertex of a path
            if (parent[i] == -1)
                root = i;
            if (path_id[i] > num_paths)
                num_paths = path_id[i];
            if (parent[i] == -1 || path[parent[i]] != i)
                for (NodeID j = i, currentPos = 0; j != -1; j = path[j])
                {
                    head[j] = i;
                    minpathPos[j] = currentPos++;
                }
        }
        num_paths++;
    }
    //function that will fill the minpath leaves with the corresponding initial weight (Karger, smallest cut that 1-respects given tree)
    // void set(NodeID v, const int &value)
    // {
    //     //TODO change s.t
    //     P.set(minpathPos[v], value);
    // }

    // int MinPath(NodeID v)
    // {
    //     int res = numeric_limits<int>::max();
    //     //go through all paths until we hit the root of spantree
    //     for (; head[v] != root; v = parent[head[v]])
    //     {
    //         res = min(res, P.minprefix(minpathPos[head[v]], minpathPos[v]));
    //     }
    //     res = min(res, P.minprefix(minpathPos[root], minpathPos[v]));
    //     return res;
    // }

    // void AddPath(NodeID v, const int &value)
    // {
    //     //go through all paths until we hit the root of spantree
    //     for (; head[v] != root; v = parent[head[v]])
    //     {
    //         P.addprefix(minpathPos[head[v]], minpathPos[v], value);
    //     }
    //     P.addprefix(minpathPos[root], minpathPos[v], value);
    // }
};

#endif