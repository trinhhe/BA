#ifndef MININUMPATH_H_
#define MININUMPATH_H_

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <string.h>
#include <assert.h>
#include <cmath>

using namespace std;

class MinimumPath
{
    //number of nodes
    int num_nodes;
    //nodes store difference of min value of right subtree and left subtree, we also store min value of complete tree in tree[0]
    pvector<int> tree;

public:
    // CONSTRUCTORS //
    MinimumPath(int n)
    {
        num_nodes = n;
        tree.reserve(2 * n);
    }

    ~MinimumPath() {}

    //function to find index in tree for some node v (range from 0 to num_nodes-1)
    // when num_nodes not power of two e.g. n=6, leaves of tree with idx 8 9 10 11 6 7, contain w1 w2 w3 w4 w5 w6 in that order
    // else e.g. n=4, leaves with idx 4 5 6 7, contain w1 w2 w3 w4 in this order
    int findIndex(int v)
    {
        assert(v < num_nodes);
        //check if num_nodes power of two
        if (num_nodes != 0 && (num_nodes & (num_nodes - 1)) == 0)
            return num_nodes + v;
        else
        {
            int left_most_node = pow(2, ceil(log2(num_nodes)));
            if (left_most_node + v < 2 * num_nodes)
                return left_most_node + v;
            else
                return left_most_node - num_nodes + v;
        }
    }

    //function to fill up leaves of MinPath
    void set(int v, int value)
    {
        tree[findIndex(v)] = value;
    }

    //initialize minpath stucture with given leave values
    void build()
    {
        //helper to store min value of each subtree
        pvector<int> store_min(2 * num_nodes);
        //store min of all leaves
        for (int i = num_nodes; i < 2 * num_nodes; ++i)
            store_min[i] = tree[i];

        for (int i = num_nodes - 1; i > 0; i--)
        {
            //min of right child - min of left child
            tree[i] = store_min[i << 1 | 1] - store_min[i << 1];
            store_min[i] = min(store_min[i << 1], store_min[i << 1 | 1]);
        }
        //store min value in first entry
        tree[0] = store_min[1];
    }

    int minprefix(NodeID node)
    {
        //d stores min weight of current subtree (where as the leaf of the weight is before the last node to be queried)
        //minus the overall min weight of current subtree
        int node_index = findIndex(node);
        int d_l = 0, d_r = 0;
        for (int i = node_index; i != 1; i >>= 1)
        {
            (i % 2 == 0) ? d_r = 0 : d_l = 0;

            //smallest weight up to queried vertex and smallest weight overall are in left subtree
            if (tree[i >> 1] > 0)
                d_r = d_l;
            //smallest weight up to queried vertex and smallest weight overall are in right subtree
            // && queried vertex must be descendant in right subtree
            else if (d_r + tree[i >> 1] < 0 && i % 2 == 1)
                d_l = d_r;
            //smallest weight up to queried vertex in left subtree and smallest weight overall in right subtree
            else
                d_r = d_l - tree[i >> 1], d_l = d_l - tree[i >> 1];
        }
        return d_l + tree[0];
    }

    void addprefix(NodeID node, const int &value)
    {
        //update all prefix nodes
        for (int i = 0; i <= node; i++)
        {
            int t = tree[findIndex(i)];
            if (value < 0)
                assert(t >= t + value);
            else
                assert(t <= t + value);
            tree[findIndex(i)] += value;
        }

        //going from last to be updated leaf to root
        //phi holds the difference between the updated min value and the previous min value of a node
        int phi_l = value, phi_r = value;
        for (int i = findIndex(node); i != 1; i >>= 1)
        {
            int delta_prev = tree[i >> 1];
            //determine missing value and compute new difference value for parent of node i
            (i % 2 == 0) ? phi_r = 0 : phi_l = value;
            tree[i >> 1] += (phi_r - phi_l);

            //determine value phi for next iteration
            if (delta_prev > 0 && tree[i >> 1] > 0)
                phi_r = phi_l;
            else if (delta_prev <= 0 && tree[i >> 1] <= 0)
                phi_l = phi_r;
            else if (delta_prev <= 0 && tree[i >> 1] > 0)
                phi_l = phi_l - delta_prev, phi_r = phi_l - delta_prev;
            else
                phi_l = phi_r + delta_prev, phi_r = phi_r + delta_prev;
        }
        //update minvalue of tree
        tree[0] += phi_l;
    }

    void print()
    {
        cout << "MinimumPath: " << endl;
        cout << "num_nodes: " << num_nodes << endl;
        for (int i = 1; i < 2 * num_nodes; i++)
        {
            cout << i << ": " << tree[i] << endl;
        }

        cout << "minvalue: " << tree[0] << endl;
    }
};

#endif