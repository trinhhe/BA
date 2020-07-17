#from https://github.com/PJK/comm-avoiding-cuts-cc

import networkx as nx
import os
from random import randint

def set_weights(G, weight):
    for (u,v,w) in G.edges(data=True):
        w['weight'] = weight
    return G

def randomize_weights(G, max):
    for (u,v,w) in G.edges(data=True):
        w['weight'] = randint(1, max)
    return G

def print_graph(G):
    # print('{0} {1}'.format(nx.number_of_nodes(G), nx.number_of_edges(G)))
    for (u,v,w) in G.edges(data=True):
        print('{0} {1} {2}'.format(u, v, w['weight']))

def read_stdin():
    first_line = input()
    while first_line.startswith('#'):
        first_line = input()
    n, m = map(int, first_line .split())
    G = nx.Graph()
    G.add_nodes_from(range(0, n))
    for _i in range(0, m):
        u, v, w = map(int, input().split())
        G.add_edge(u, v, weight=w)
    return G
