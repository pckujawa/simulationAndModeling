#-------------------------------------------------------------------------------
# Name:        graph library
# Author:      Pat Kujawa
#-------------------------------------------------------------------------------
#!/usr/bin/env python
##from __future__ import division
##from scipy.integrate import odeint  # for integrate.odeint
##import numpy as np
##import pylab as pl
##import math
import unittest
##import itertools
##from pprint import pprint, pformat
from collections import defaultdict

##from libs import Struct

class GraphList(object):
    def __init__(self, edges='undirected'):
        """
        """
        self._edges = edges
        self._adj = defaultdict(list)

    def add(self, index, *neighbor_ixs):
        self._adj[index].extend(neighbor_ixs)
        if self._edges == 'undirected':
            for n in neighbor_ixs:
                self._adj[n].append(index)
        return self

    def neighbors(self, index):
        """Return a copy of the neighbors of the specified item.
        """
        ns = self._adj[index][:]
        return ns


class GraphListUndirectedTests(unittest.TestCase):
    def get_target(self):
        return GraphList(edges='undirected')

    def test_neighbor_of_added_item_ix_is_second_added_item(self):
        target = self.get_target()
        target.add(0, 1)
        expected = [0]
        actual = target.neighbors(1)
        self.assertEqual(expected, actual)

    def test_given_undirected_edge_then_connected_nodes_are_neighbors(self):
        target = self.get_target()
        target.add(1, 2)
        expected = {1: [2], 2: [1]}
        actual = {1: target.neighbors(1), 2: target.neighbors(2)}
        self.assertEqual(expected, actual)


class GraphListDirectedTests(unittest.TestCase):
    def get_target(self):
        return GraphList(edges='directed')

    def test_given_directed_edge_then_only_endpoint_node_is_neighbor(self):
        target = self.get_target()
        target.add(1, 2)
        expected = {1: [2], 2: []}
        actual = {1: target.neighbors(1), 2: target.neighbors(2)}
        self.assertEqual(expected, actual)


if __name__ == '__main__':
    unittest.main(verbosity=2)
