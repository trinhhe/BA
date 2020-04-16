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
    }

    //initialize minpath stucture with given leave values
    void build()
    {
    }

    int minprefix(NodeID v)
    {
        return 0;
    }

    void addprefix(NodeID v, const int &value)
    {
    }
};

#endif