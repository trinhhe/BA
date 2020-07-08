#ifndef TREESGENERATOR_H_
#define TREESGENERATOR_H_

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

#include "graph.h"
#include "pvector.h"
#include "builder.h"
#include "util.h"
#include "timer.h"

/*
Author: Henry Trinh

class to represent treesgenerator
*/
using namespace std;
typedef int32_t NodeID;
typedef int32_t WeightT;
typedef NodeWeight<NodeID, WeightT> WNode;
typedef EdgePair<NodeID, WNode> WEdge;

const double max_allowed_deviation = 1e-10;

class TreesGenerator
{
    const pvector<WEdge> &G;

    double d, eps1, eps2, f, b;
    int n, m, weight_cap, c_dash;
};

#endif