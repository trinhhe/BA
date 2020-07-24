#!/usr/bin/env python3
#from https://github.com/PJK/comm-avoiding-cuts-cc

import networkx as nx
import sys, os
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from graph_generator import io_utils
import argparse
import random, datetime

parser = argparse.ArgumentParser(
	description='Generate a graph. For available generators and usage, see https://networkx.github.io/documentation/networkx-1.10/reference/generators.html',
	formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument('-w', '--weight', metavar='W', type=int, help='edge weight', default=100)
parser.add_argument('-s', '--seed', metavar='S', type=int, help='randomization seed', default=1234)
parser.add_argument('--randomize', dest='set_weights', action='store_const',
					const=io_utils.randomize_weights, default=io_utils.set_weights,
					help='randomize edge weights uniformly within [1, W)')
parser.add_argument('--randomize_normal', dest='set_weights', action='store_const',
                    const=io_utils.randomize_weights_normal_distribution, default=io_utils.set_weights,
                    help='randomize edge weights normally with W as mean')
parser.add_argument('graph', metavar='G', type=str, help='graph specification (e.g. complete_graph(N), ...)')
parser.add_argument('N', metavar='N', type=int, help='vertex count')
parser.add_argument('-m', '--edges', metavar='M', type=int, help='edge count', nargs='?')
parser.add_argument('-p', '--prob', metavar='P', type=float, help='edge creation probability (graph-specific)', nargs='?')
parser.add_argument('-k', '--k', metavar='K', type=int, help='neighbor distance/connectivity/clusters count (graph-specific)', nargs='?')
parser.add_argument('-nc', '--num_cliques', metavar='NC', type=int, help='num_cliques (graph specific for ring_of_cliques)', nargs='?')
parser.add_argument('-cs', '--clique_size', metavar='CS', type=int, help='clique_size (graph specific for ring_of_cliques)', nargs='?')
parser.add_argument('-pl', '--partition_list', metavar='PL', type=int, help='size of groups (graph specific for random_partition',nargs='*')
parser.add_argument('-pin', '--p_in', metavar='PIN', type=float, help='probability of edges with in groups (graph specific for random_partition',nargs='?')
parser.add_argument('-pout', '--p_out', metavar='POUT', type=float, help='probability of edges between groups (graph specific for random_partition',nargs='?')

args = parser.parse_args()
random.seed(args.seed)

G = eval('nx.' + args.graph, {**globals(), **{'N': args.N, 'M': args.edges, 'P': args.prob, 'K': args.k, 'NC': args.num_cliques,'CS': args.clique_size, 'PL': args.partition_list, 'PIN': args.p_in, 'POUT': args.p_out}}, locals())

# print('# {} {} {}'.format(datetime.datetime.now(), os.popen('git rev-parse HEAD').read().strip(), args))
args.set_weights(G, args.weight)
io_utils.print_graph(G)
