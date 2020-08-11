#include <algorithm>
#include <iostream>
#include <assert.h>


#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"
#include "timer.h"
#include "util.h"
#include "spantree.h"
#include "minimumpath.h"

int testComparablecut1()
{
    pvector<WEdge> g(9);
    g[0] = WEdge(0,WNode(1,1));
    g[1] = WEdge(0,WNode(2,1));
    g[2] = WEdge(0,WNode(4,1));
    g[3] = WEdge(1,WNode(2,1));
    g[4] = WEdge(1,WNode(3,1));
    g[5] = WEdge(1,WNode(4,1));
    g[6] = WEdge(2,WNode(5,1));
    g[7] = WEdge(3,WNode(5,1));
    g[8] = WEdge(4,WNode(5,1));
    pvector<WEdge> t(5);
    t[0] = WEdge(0,WNode(1,1));
    t[1] = WEdge(1,WNode(2,1));
    t[2] = WEdge(1,WNode(3,1));
    t[3] = WEdge(1,WNode(4,1));
    t[4] = WEdge(3,WNode(5,1));
    auto G =  WeightedBuilder::Load_CSR_From_Edgelist(g, true);
    auto T =  WeightedBuilder::Load_CSR_From_Edgelist(t, true);
    int root = 0;
    int num_nodes = 6;
    int num_paths = 4;
    pvector<NodeID> parent(6);
    parent[0] = -1;
    parent[1] = 0;
    parent[2] = 1;
    parent[3] = 1;
    parent[4] = 1;
    parent[5] = 3;
    
    pvector<NodeID> path_id(6);
    path_id[0] = 0;
    path_id[1] = 0;
    path_id[2] = 1;
    path_id[3] = 2;
    path_id[4] = 3;
    path_id[5] = 2;

    pvector<NodeID> path(6);
    path[0] = 1;
    path[1] = -1;
    path[2] = -1;
    path[3] = 5;
    path[4] = -1;
    path[5] = -1;

    pvector<NodeID> head(6);
    head[0] = 0;
    head[1] = 0;
    head[2] = 2;
    head[3] = 3;
    head[4] = 4;
    head[5] = 3;

    pvector<NodeID> minpathPos(6);
    minpathPos[0] = 0;
    minpathPos[1] = 1;
    minpathPos[2] = 0;
    minpathPos[3] = 0;
    minpathPos[4] = 0;
    minpathPos[5] = 1;

    pvector<PathID> lengths(4);
    lengths[0] = 2;
    lengths[1] = 1;
    lengths[2] = 2;
    lengths[3] = 1;

    SpanTree *st = new SpanTree(G,T,root,num_nodes,num_paths,parent, path_id,path,head,minpathPos, lengths);
    int which_sol;
    return st->compute(which_sol);  
}

int testComparablecut2()
{
    pvector<WEdge> g(5);
    g[0] = WEdge(0,WNode(1,1));
    g[1] = WEdge(0,WNode(4,2));
    g[2] = WEdge(1,WNode(2,2));
    g[3] = WEdge(2,WNode(3,1));
    g[4] = WEdge(3,WNode(4,3));
    pvector<WEdge> t(4);
    t[0] = WEdge(0,WNode(1,1));
    t[1] = WEdge(1,WNode(2,2));
    t[2] = WEdge(2,WNode(3,1));
    t[3] = WEdge(3,WNode(4,3));
    auto G =  WeightedBuilder::Load_CSR_From_Edgelist(g, true);
    auto T =  WeightedBuilder::Load_CSR_From_Edgelist(t, true);
    int root = 0;
    int num_nodes = 5;
    int num_paths = 1;
    pvector<NodeID> parent(5);
    parent[0] = -1;
    parent[1] = 0;
    parent[2] = 1;
    parent[3] = 2;
    parent[4] = 3;
    
    pvector<NodeID> path_id(5);
    path_id[0] = 0;
    path_id[1] = 0;
    path_id[2] = 0;
    path_id[3] = 0;
    path_id[4] = 0;

    pvector<NodeID> path(5);
    path[0] = 1;
    path[1] = 2;
    path[2] = 3;
    path[3] = 4;
    path[4] = -1;

    pvector<NodeID> head(5);
    head[0] = 0;
    head[1] = 0;
    head[2] = 0;
    head[3] = 0;
    head[4] = 0;

    pvector<NodeID> minpathPos(5);
    minpathPos[0] = 0;
    minpathPos[1] = 1;
    minpathPos[2] = 2;
    minpathPos[3] = 3;
    minpathPos[4] = 4;

    pvector<PathID> lengths(1);
    lengths[0] = 5;
    

    SpanTree *st = new SpanTree(G,T,root,num_nodes,num_paths,parent, path_id,path,head,minpathPos, lengths);
    int which_sol;
    return st->compute(which_sol);  
}

int testComparablecut3()
{
    pvector<WEdge> g(8);
    g[0] = WEdge(0,WNode(5,4));
    g[1] = WEdge(0,WNode(6,5));
    g[2] = WEdge(0,WNode(1,1));
    g[3] = WEdge(1,WNode(3,5));
    g[4] = WEdge(1,WNode(2,6));
    g[5] = WEdge(2,WNode(4,1));
    g[6] = WEdge(3,WNode(4,1));
    g[7] = WEdge(4,WNode(6,7));
    pvector<WEdge> t(6);
    t[0] = WEdge(0,WNode(5,4));
    t[1] = WEdge(0,WNode(6,5));
    t[2] = WEdge(0,WNode(1,1));
    t[3] = WEdge(1,WNode(3,5));
    t[4] = WEdge(1,WNode(2,6));
    t[5] = WEdge(2,WNode(4,1));

    auto G =  WeightedBuilder::Load_CSR_From_Edgelist(g, true);
    auto T =  WeightedBuilder::Load_CSR_From_Edgelist(t, true);
    int root = 0;
    int num_nodes = 7;
    int num_paths = 5;
    pvector<NodeID> parent(7);
    parent[0] = -1;
    parent[1] = 0;
    parent[2] = 1;
    parent[3] = 1;
    parent[4] = 2;
    parent[5] = 0;
    parent[6] = 0;
    
    pvector<NodeID> path_id(7);
    path_id[0] = 0;
    path_id[1] = 0;
    path_id[2] = 1;
    path_id[3] = 2;
    path_id[4] = 1;
    path_id[5] = 3;
    path_id[6] = 4;

    pvector<NodeID> path(7);
    path[0] = 1;
    path[1] = -1;
    path[2] = 4;
    path[3] = -1;
    path[4] = -1;
    path[5] = -1;
    path[6] = -1;

    pvector<NodeID> head(7);
    head[0] = 0;
    head[1] = 0;
    head[2] = 2;
    head[3] = 3;
    head[4] = 2;
    head[5] = 5;
    head[6] = 6;

    pvector<NodeID> minpathPos(7);
    minpathPos[0] = 0;
    minpathPos[1] = 1;
    minpathPos[2] = 0;
    minpathPos[3] = 0;
    minpathPos[4] = 1;
    minpathPos[5] = 0;
    minpathPos[6] = 0;

    pvector<PathID> lengths(5);
    lengths[0] = 2;
    lengths[1] = 2;
    lengths[2] = 1;
    lengths[3] = 1;
    lengths[4] = 1;

    SpanTree *st = new SpanTree(G,T,root,num_nodes,num_paths,parent, path_id,path,head,minpathPos, lengths);
    int which_sol;
    return st->compute(which_sol);  
}

int testComparablecut4()
{
    pvector<WEdge> g(9);
    g[0] = WEdge(0,WNode(3,8));
    g[1] = WEdge(0,WNode(4,7));
    g[2] = WEdge(0,WNode(1,1));
    g[3] = WEdge(1,WNode(2,1));
    g[4] = WEdge(4,WNode(8,3));
    g[5] = WEdge(2,WNode(5,4));
    g[6] = WEdge(2,WNode(6,3));
    g[7] = WEdge(6,WNode(7,4));
    g[8] = WEdge(6,WNode(8,5));
    pvector<WEdge> t(8);
    t[0] = WEdge(0,WNode(3,8));
    t[1] = WEdge(0,WNode(4,7));
    t[2] = WEdge(0,WNode(1,1));
    t[3] = WEdge(1,WNode(2,1));
    t[4] = WEdge(2,WNode(5,4));
    t[5] = WEdge(2,WNode(6,3));
    t[6] = WEdge(6,WNode(7,4));
    t[7] = WEdge(6,WNode(8,5));
    auto G =  WeightedBuilder::Load_CSR_From_Edgelist(g, true);
    auto T =  WeightedBuilder::Load_CSR_From_Edgelist(t, true);
    int root = 0;
    int num_nodes = 9;
    int num_paths = 6;
    pvector<NodeID> parent(9);
    parent[0] = -1;
    parent[1] = 0;
    parent[2] = 1;
    parent[3] = 0;
    parent[4] = 0;
    parent[5] = 2;
    parent[6] = 2;
    parent[7] = 6;
    parent[8] = 6;
    
    pvector<NodeID> path_id(9);
    path_id[0] = 0;
    path_id[1] = 0;
    path_id[2] = 0;
    path_id[3] = 1;
    path_id[4] = 2;
    path_id[5] = 3;
    path_id[6] = 0;
    path_id[7] = 4;
    path_id[8] = 5;

    pvector<NodeID> path(9);
    path[0] = 1;
    path[1] = 2;
    path[2] = 6;
    path[3] = -1;
    path[4] = -1;
    path[5] = -1;
    path[6] = -1;
    path[7] = -1;
    path[8] = -1;
    pvector<NodeID> head(9);
    head[0] = 0;
    head[1] = 0;
    head[2] = 0;
    head[3] = 3;
    head[4] = 4;
    head[5] = 5;
    head[6] = 0;
    head[7] = 7;
    head[8] = 8;

    pvector<NodeID> minpathPos(9);
    minpathPos[0] = 0;
    minpathPos[1] = 1;
    minpathPos[2] = 2;
    minpathPos[3] = 0;
    minpathPos[4] = 0;
    minpathPos[5] = 0;
    minpathPos[6] = 3;
    minpathPos[7] = 0;
    minpathPos[8] = 0;

    pvector<PathID> lengths(9);
    lengths[0] = 4;
    lengths[1] = 1;
    lengths[2] = 1;
    lengths[3] = 1;
    lengths[4] = 1;
    lengths[5] = 1;
    

    SpanTree *st = new SpanTree(G,T,root,num_nodes,num_paths,parent, path_id,path,head,minpathPos, lengths);
    int which_sol;
    return st->compute(which_sol);  
}

int testComparablecut5()
{
    pvector<WEdge> g(13);
    g[0] = WEdge(0,WNode(3,8));
    g[1] = WEdge(0,WNode(4,7));
    g[2] = WEdge(0,WNode(1,1));
    g[3] = WEdge(1,WNode(2,1));
    g[4] = WEdge(3,WNode(5,1));
    g[5] = WEdge(4,WNode(8,1));
    g[6] = WEdge(2,WNode(5,4));
    g[7] = WEdge(2,WNode(6,3));
    g[8] = WEdge(5,WNode(7,3));
    g[9] = WEdge(6,WNode(7,4));
    g[10] = WEdge(6,WNode(8,5));
    g[11] = WEdge(5,WNode(9,5));
    g[12] = WEdge(5,WNode(10,5));
    pvector<WEdge> t(10);
    t[0] = WEdge(0,WNode(3,8));
    t[1] = WEdge(0,WNode(4,7));
    t[2] = WEdge(0,WNode(1,1));
    t[3] = WEdge(1,WNode(2,1));
    t[4] = WEdge(2,WNode(5,4));
    t[5] = WEdge(2,WNode(6,3));
    t[6] = WEdge(6,WNode(7,4));
    t[7] = WEdge(6,WNode(8,5));
    t[8] = WEdge(5, WNode(9,5));
    t[9] = WEdge(5, WNode(10,5));
    auto G =  WeightedBuilder::Load_CSR_From_Edgelist(g, true);
    auto T =  WeightedBuilder::Load_CSR_From_Edgelist(t, true);
    int root = 0;
    int num_nodes = 11;
    int num_paths = 9;
    pvector<NodeID> parent(11);
    parent[0] = -1;
    parent[1] = 0;
    parent[2] = 1;
    parent[3] = 0;
    parent[4] = 0;
    parent[5] = 2;
    parent[6] = 2;
    parent[7] = 6;
    parent[8] = 6;
    parent[9] = 5;
    parent[10] = 5;
    
    pvector<NodeID> path_id(11);
    path_id[0] = 0;
    path_id[1] = 0;
    path_id[2] = 0;
    path_id[3] = 1;
    path_id[4] = 2;
    path_id[5] = 3;
    path_id[6] = 4;
    path_id[7] = 5;
    path_id[8] = 6;
    path_id[9] = 7;
    path_id[10] = 8;

    pvector<NodeID> path(11);
    path[0] = 1;
    path[1] = 2;
    path[2] = -1;
    path[3] = -1;
    path[4] = -1;
    path[5] = -1;
    path[6] = -1;
    path[7] = -1;
    path[8] = -1;
    path[9] = -1;
    path[10] = -1;
    pvector<NodeID> head(11);
    head[0] = 0;
    head[1] = 0;
    head[2] = 0;
    head[3] = 3;
    head[4] = 4;
    head[5] = 5;
    head[6] = 6;
    head[7] = 7;
    head[8] = 8;
    head[9] = 9;
    head[10] = 10;

    pvector<NodeID> minpathPos(11);
    minpathPos[0] = 0;
    minpathPos[1] = 1;
    minpathPos[2] = 2;
    minpathPos[3] = 0;
    minpathPos[4] = 0;
    minpathPos[5] = 0;
    minpathPos[6] = 0;
    minpathPos[7] = 0;
    minpathPos[8] = 0;
    minpathPos[9] = 0;
    minpathPos[10] = 0;

    pvector<PathID> lengths(9);
    lengths[0] = 3;
    lengths[1] = 1;
    lengths[2] = 1;
    lengths[3] = 1;
    lengths[4] = 1;
    lengths[5] = 1;
    lengths[6] = 1; 
    lengths[7] = 1;
    lengths[8] = 1;

    SpanTree *st = new SpanTree(G,T,root,num_nodes,num_paths,parent, path_id,path,head,minpathPos, lengths);
    int which_sol;
    return st->compute(which_sol);  
}

int testIncomparablecut1()
{
    pvector<WEdge> g(4);
    g[0] = WEdge(0,WNode(1,4));
    g[1] = WEdge(1,WNode(2,1));
    g[2] = WEdge(1,WNode(3,1));
    g[3] = WEdge(2,WNode(3,3));
    
    pvector<WEdge> t(3);
    t[0] = WEdge(0,WNode(1,4));
    t[1] = WEdge(1,WNode(2,1));
    t[2] = WEdge(1,WNode(3,1));

    auto G =  WeightedBuilder::Load_CSR_From_Edgelist(g, true);
    auto T =  WeightedBuilder::Load_CSR_From_Edgelist(t, true);
    int root = 0;
    int num_nodes = 4;
    int num_paths = 3;
    pvector<NodeID> parent(4);
    parent[0] = -1;
    parent[1] = 0;
    parent[2] = 1;
    parent[3] = 1;
    pvector<NodeID> path_id(4);
    path_id[0] = 0;
    path_id[1] = 0;
    path_id[2] = 1;
    path_id[3] = 2;
    pvector<NodeID> path(4);
    path[0] = 1;
    path[1] = -1;
    path[2] = -1;
    path[3] = -1;
    pvector<NodeID> head(4);
    head[0] = 0;
    head[1] = 0;
    head[2] = 2;
    head[3] = 3;
    pvector<NodeID> minpathPos(4);
    minpathPos[0] = 0;
    minpathPos[1] = 1;
    minpathPos[2] = 0;
    minpathPos[3] = 0;
    pvector<PathID> lengths(3);
    lengths[0] = 2;
    lengths[1] = 1;
    lengths[2] = 1;

    SpanTree *st = new SpanTree(G,T,root,num_nodes,num_paths,parent, path_id,path,head,minpathPos, lengths);
    int which_sol;
    return st->compute(which_sol);  
}

int testIncomparablecut2()
{
    pvector<WEdge> g(7);
    g[0] = WEdge(0,WNode(1,5));
    g[1] = WEdge(1,WNode(3,4));
    g[2] = WEdge(1,WNode(4,1));
    g[3] = WEdge(1,WNode(5,4));
    g[4] = WEdge(0,WNode(2,1));
    g[5] = WEdge(2,WNode(5,1));
    g[6] = WEdge(2,WNode(4,5));
    
    pvector<WEdge> t(5);
    t[0] = WEdge(0,WNode(1,5));
    t[1] = WEdge(1,WNode(3,4));
    t[2] = WEdge(1,WNode(4,1));
    t[3] = WEdge(1,WNode(5,4));
    t[4] = WEdge(0,WNode(2,1));

    auto G =  WeightedBuilder::Load_CSR_From_Edgelist(g, true);
    auto T =  WeightedBuilder::Load_CSR_From_Edgelist(t, true);
    int root = 0;
    int num_nodes = 6;
    int num_paths = 5;
    pvector<NodeID> parent(6);
    parent[0] = -1;
    parent[1] = 0;
    parent[2] = 0;
    parent[3] = 1;
    parent[4] = 1;
    parent[5] = 1;
    pvector<NodeID> path_id(6);
    path_id[0] = 0;
    path_id[1] = 0;
    path_id[2] = 1;
    path_id[3] = 2;
    path_id[4] = 3;
    path_id[5] = 4;
    pvector<NodeID> path(6);
    path[0] = 1;
    path[1] = -1;
    path[2] = -1;
    path[3] = -1;
    path[4] = -1;
    path[5] = -1;
    pvector<NodeID> head(6);
    head[0] = 0;
    head[1] = 0;
    head[2] = 2;
    head[3] = 3;
    head[4] = 4;
    head[5] = 5;
    pvector<NodeID> minpathPos(6);
    minpathPos[0] = 0;
    minpathPos[1] = 1;
    minpathPos[2] = 0;
    minpathPos[3] = 0;
    minpathPos[4] = 0;
    minpathPos[5] = 0;
    pvector<PathID> lengths(5);
    lengths[0] = 2;
    lengths[1] = 1;
    lengths[2] = 1;
    lengths[3] = 1;
    lengths[4] = 1;

    SpanTree *st = new SpanTree(G,T,root,num_nodes,num_paths,parent, path_id,path,head,minpathPos, lengths);
    int which_sol;
    return st->compute(which_sol);  
}

int testIncomparablecut3()
{
    pvector<WEdge> g(10);
    g[0] = WEdge(0,WNode(1,1));
    g[1] = WEdge(0,WNode(2,1));
    g[2] = WEdge(1,WNode(3,3));
    g[3] = WEdge(1,WNode(4,3));
    g[4] = WEdge(2,WNode(5,3));
    g[5] = WEdge(2,WNode(6,3));
    g[6] = WEdge(2,WNode(7,3));
    g[7] = WEdge(2,WNode(8,3));
    g[8] = WEdge(1,WNode(2,4));
    g[9] = WEdge(4,WNode(5,6));

    pvector<WEdge> t(8);
    t[0] = WEdge(0,WNode(1,1));
    t[1] = WEdge(0,WNode(2,1));
    t[2] = WEdge(1,WNode(3,3));
    t[3] = WEdge(1,WNode(4,3));
    t[4] = WEdge(2,WNode(5,3));
    t[5] = WEdge(2,WNode(6,3));
    t[6] = WEdge(2,WNode(7,3));
    t[7] = WEdge(2,WNode(8,3));

    auto G =  WeightedBuilder::Load_CSR_From_Edgelist(g, true);
    auto T =  WeightedBuilder::Load_CSR_From_Edgelist(t, true);
    int root = 0;
    int num_nodes = 9;
    int num_paths = 9;
    pvector<NodeID> parent(9);
    parent[0] = -1;
    parent[1] = 0;
    parent[2] = 0;
    parent[3] = 1;
    parent[4] = 1;
    parent[5] = 2;
    parent[6] = 2;
    parent[7] = 2;
    parent[8] = 2;

    pvector<NodeID> path_id(9);
    path_id[0] = 0;
    path_id[1] = 1;
    path_id[2] = 2;
    path_id[3] = 3;
    path_id[4] = 4;
    path_id[5] = 5;
    path_id[6] = 6;
    path_id[7] = 7;
    path_id[8] = 8;
    path_id[9] = 9;

    pvector<NodeID> path(9);
    path[0] = -1;
    path[1] = -1;
    path[2] = -1;
    path[3] = -1;
    path[4] = -1;
    path[5] = -1;
    path[6] = -1;
    path[7] = -1;
    path[8] = -1;

    pvector<NodeID> head(9);
    head[0] = 0;
    head[1] = 1;
    head[2] = 2;
    head[3] = 3;
    head[4] = 4;
    head[5] = 5;
    head[6] = 6;
    head[7] = 7;
    head[8] = 8;

    pvector<NodeID> minpathPos(9);
    minpathPos[0] = 0;
    minpathPos[1] = 0;
    minpathPos[2] = 0;
    minpathPos[3] = 0;
    minpathPos[4] = 0;
    minpathPos[5] = 0;
    minpathPos[6] = 0;
    minpathPos[7] = 0;
    minpathPos[8] = 0;

    pvector<PathID> lengths(9);
    lengths[0] = 1;
    lengths[1] = 1;
    lengths[2] = 1;
    lengths[3] = 1;
    lengths[4] = 1;
    lengths[5] = 1;
    lengths[6] = 1;
    lengths[7] = 1;
    lengths[8] = 1;
    SpanTree *st = new SpanTree(G,T,root,num_nodes,num_paths,parent, path_id,path,head,minpathPos, lengths);
    int which_sol;
    return st->compute(which_sol);  
}
int main(int argc, char *argv[])
{
    assert(testComparablecut1() == 2);
    assert(testComparablecut2() == 2);
    assert(testComparablecut3() == 3);
    assert(testComparablecut4() == 2);
    assert(testComparablecut5() == 2);
    assert(testIncomparablecut1()== 2);
    assert(testIncomparablecut2()==3);
    assert(testIncomparablecut3()==2);

    return 0;

}