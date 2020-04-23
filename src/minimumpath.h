#ifndef MININUMPATH_H_
#define MININUMPATH_H_

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <string.h>
#include <assert.h>

class MinimumPath
{
    //number of nodes
    int num_nodes;
    pvector<int> tree;

public:
    // CONSTRUCTORS //
    MinimumPath(int n)
    {
        num_nodes = n;
        tree.reserve(2 * n);
    }

    ~MinimumPath() {}

    //function to fill up leaves of MinPath
    void set(int v, const int &value)
    {
        tree[num_nodes + v] = value;
    }

    //initialize minpath stucture with given leave values
    void build()
    {
        for (int i = num_nodes - 1; i > 0; i--)
        {
                }
    }

    int minprefix(NodeID startnode, NodeID lastnode)
    {
        return 0;
    }

    void addprefix(NodeID startnode, NodeID lastnode, const int &value)
    {
    }
};

#endif