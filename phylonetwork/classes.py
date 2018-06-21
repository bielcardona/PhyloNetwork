from networkx import DiGraph, is_directed_acyclic_graph, dfs_successors
from networkx import single_source_shortest_path_length, all_pairs_shortest_path_length, dijkstra_path_length
# from networkx import pydot_layout, draw_networkx
import networkx as nx

import numpy, pyparsing, copy

from .eNewick import eNewickParser
from .utils import total_cmp, cmp_to_key
from . import permutations
from .memoize import memoize_method
# def memoize_method(f): return f # use it to document memoized methods (sphinx bug?)
from .exceptions import MalformedNewickException

from itertools import combinations
from collections import Counter

import logging
logger = logging.getLogger(__name__)


class PhyloNetwork(DiGraph):
    """
    Main class for phylogenetic networks and trees (with or without nested taxa).
    
    You can create a PhyloNetwork with its eNewick representation::
    
        >>> network = PhyloNetwork(eNewick="((1,2),(3,4)5);")
        >>> network.taxa()
        ... ['1','2','3','4','5']
        >>> network.leaves()
        ... ['_3', '_4', '_6', '_7']
        >>> network.label('_6')
        ... '3'
      
    If your eNewick string is malformed, you'll receive a MalformedNewickException::
    
        >>> network = PhyloNetwork(eNewick="(1)")
        ... Traceback (most recent call last):
        ... (...)
        ... PhyloNetwork.classes.MalformedNewickException: Malformed eNewick string
        
    You can also start with an existing networkx graph::
    
        >>> graph = networkx.DiGraph()
        >>> graph.add_nodes_from(range(5))
        >>> graph.add_edges_from([(0,1), (0,2), (1,3), (1,4)])
        >>> network = PhyloNetwork(data = graph)
        >>> network.leaves()
        ... [2, 3, 4]
        >>> network.taxa()
        ... []
        
    """

    def __init__(self, data=None, name='', eNewick=None, ignore_prefix=None, id_offset=0):
        # initialization here
        DiGraph.__init__(self, data)
        self.name = name
        self._labels = {}
        self._lastlabel = id_offset
        self.cache = {}
        if eNewick != None:
            self._from_eNewick(eNewick, ignore_prefix=ignore_prefix)

    @memoize_method
    def is_phylogenetic_network(self):
        """
        Returns True if the network is a Phylogenetic Network. False otherwise.
                
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2),(3,4)5);")
            >>> network.is_phylogenetic_network()
            ... True
            >>> graph = networkx.DiGraph()
            >>> graph.add_nodes_from(range(4))
            >>> graph.add_edges_from([(0,1), (0,2), (1,2), (2,3), (3,1)])
            >>> PhyloNetwork(data=graph).is_phylogenetic_network()
            ... False
            
        """
        if not is_directed_acyclic_graph(self):
            return False
        return True

    def set_label(self, node, label):
        """
        Set a new label to a node.
        
        EXAMPLE::
    
            >>> graph = networkx.DiGraph()
            >>> graph.add_nodes_from(range(5))
            >>> graph.add_edges_from([(0,1), (0,2), (1,3), (1,4)])
            >>> network = PhyloNetwork(data = graph)
            >>> network.leaves()
            ... [2, 3, 4]
            >>> network.taxa()
            ... []
            >>> network.set_label(2, 'Label 1')
            >>> network.set_label(3, 'Label 2')
            >>> network.set_label(4, 'Label 3')
            >>> network.taxa()
            ... ['Label 1', 'Label 2', 'Label 3']
            >>> network.set_label(2, 'New label')
            >>> network.taxa()
            ... ['Label 2', 'Label 3', 'New label']

        """
        if node in self:
            self._labels[node] = label
            self.cache = {}

    @memoize_method
    def taxa(self):
        """
        Returns the taxa (set of labels) of self.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2),(3,4)5);")
            >>> network.taxa()
            ... ['1', '2', '3', '4', '5']
            >>> network.leaves()
            ... ['_3', '_4', '_6', '_7']
            >>> network.label('_6')
            ... '3'
            >>> network.nodes()
            ... ['_7', '_6', '_5', '_4', '_3', '_2', '_1']
            
        """
        taxa = list(set(self._labels.values()))
        taxa.sort()
        return taxa

    def label(self, node):
        """
        Returns the label of node, or None if not labelled.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2),(3,4)5);")
            >>> network.taxa()
            ... ['1', '2', '3', '4', '5']
            >>> network.leaves()
            ... ['_3', '_4', '_6', '_7']
            >>> network.label('_6')
            ... '3'
            >>> network.nodes()
            ... ['_7', '_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.label('_1')
            ... None
            
        """
        return self._labels.get(node)

    @memoize_method
    def node_by_taxa(self, taxa):
        """
        Returns the node labelled by taxa or None if no node is labelled by taxa.
        Important: If more than one node is labelled with taxa, only the first one will be displayed. 
        In order to get all nodes with a fixed label, use nodes_by_taxa.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2),(3,4)1);")
            >>> network.taxa()
            ... ['1', '2', '3', '4']
            >>> network.nodes()
            ... ['_7', '_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.node_by_taxa('3')
            ... '_6'
            >>> network.node_by_taxa('non-existing taxa')
            ... None
            >>> network.node_by_taxa('1')
            ... '_5'
            
        """
        for node in self.labelled_nodes():
            if self.label(node) == taxa:
                return node
        return None

    @memoize_method
    def nodes_by_taxa(self, taxa):
        """
        Returns all nodes labelled with taxa.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2),(3,4)1);")
            >>> network.taxa()
            ... ['1', '2', '3', '4']
            >>> network.nodes()
            ... ['_7', '_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.node_by_taxa('3')
            ... '_6'
            >>> network.nodes_by_taxa('3')
            ... ['_6']
            >>> network.node_by_taxa('1')
            ... '_5'
            >>> network.nodes_by_taxa('1')
            ... ['_5', '_3']
            >>> network.nodes_by_taxa('non-existing taxa')
            ... set([])
            
        """
        tmp = []
        for node in self.labelled_nodes():
            if self.label(node) == taxa:
                tmp.append(node)
        return set(tmp)

    def is_tree_node(self, u):
        """
        Returns True if u is a tree node, False otherwise.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((3),(1,2));")
            >>> network.nodes()
            ... ['_6', '_5', '_4', '_3', '_2', '_1']
            >>> list(map(network.is_tree_node, network.nodes()))
            ... [True, True, True, True, True, True]
            >>> network.is_tree_node('non-existing node')
            ... False
            >>> network = PhyloNetwork(eNewick="((4,5#1)2,(#1,6)3)1;")
            >>> network.nodes()
            ... ['_5', '_4', '_3', '_2', '_1', '#1']
            >>> list(map(network.is_tree_node, network.nodes()))
            ... [True, True, True, True, True, False]
            
        """
        return u in self.nodes() and self.in_degree(u) <= 1

    def is_hybrid_node(self, u):
        """
        Returns True if u is not a tree node, False otherwise.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((3),(1,2));")
            >>> network.nodes()
            ... ['_6', '_5', '_4', '_3', '_2', '_1']
            >>> list(map(network.is_hybrid_node, network.nodes()))
            ... [False, False, False, False, False, False]
            >>> network.is_hybrid_node('non-existing node')
            ... False
            >>> network = PhyloNetwork(eNewick="((4,5#1)2,(#1,6)3)1;")
            >>> network.nodes()
            ... ['_5', '_4', '_3', '_2', '_1', '#1']
            >>> list(map(network.is_hybrid_node, network.nodes()))
            ... [False, False, False, False, False, True]
            
        """
        return u in self.nodes() and self.in_degree(u) > 1

    def is_leaf(self, u):
        """
        Returns True if u is a leaf, False otherwise.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((3),(1,2));")
            >>> network.nodes()
            ... ['_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.leaves()
            ... ['_3', '_5', '_6']
            >>> network.is_leaf('_1')
            ... False
            >>> network.is_leaf('_6')
            ... True
            >>> network.is_leaf('non-existing node')
            ... False
            
        """
        return u in self.nodes() and self.out_degree(u) == 0

    def is_root(self, u):
        """
        Returns True if u is a root, False otherwise.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((3),(1,2));")
            >>> network.nodes()
            ... ['_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.leaves()
            ... ['_3', '_5', '_6']
            >>> network.is_root('_3')
            ... False
            >>> network.is_root('_1')
            ... True
            
        """
        return u in self.nodes() and self.in_degree(u) == 0

    def is_elementary_node(self, u):
        """
        Returns True if u is an elementary node, False otherwise.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1)E,2);")
            >>> network.nodes()
            ... ['_4', '_3', '_2', '_1']
            >>> map(network.is_elementary_node, network.nodes())
            ... [False, False, True, False]
            >>> network.label('_2')
            ... 'E'
            >>> network.elementary_nodes()
            ... ['_2']
            
        """

        return u in self.nodes() and ((self.in_degree(u) <= 1) and (self.out_degree(u) == 1))

    def is_labelled(self, u):
        """
        Returns True if u is a labelled node, False otherwise.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((3),(1,2))4;")
            >>> network.nodes()
            ... ['_6', '_5', '_4', '_3', '_2', '_1']
            >>> list(map(network.is_labelled, network.nodes()))
            ... [True, True, False, True, False, True]
            >>> network.label('_1')
            ... '4'
        
        """
        return u in self._labels

    @memoize_method
    def leaves(self):
        """
        Returns the set of leaves of self.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((3),(1,(5,6)2))4;")
            >>> network.nodes()
            ... ['_8', '_7', '_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.leaves()
            ... ['_3', '_5', '_7', '_8']
            >>> list(map(network.is_leaf, network.leaves()))
            ... [True, True, True, True]
            
        """
        leaves = list(filter(self.is_leaf, self.nodes()))
        leaves.sort()
        return leaves

    @memoize_method
    def roots(self):
        """
        Returns the set of roots of self.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="(1,2,3)ROOT;")
            >>> network.nodes()
            ... ['_4', '_3', '_2', '_1']
            >>> network.roots()
            ... ['_1']
            >>> network.label('_1')
            ... 'ROOT'
            
        EXAMPLE::
        
        >>> graph = networkx.DiGraph()
        >>> graph.add_nodes_from(range(5))
        >>> graph.add_edges_from([(0,1), (2,3), (1,4), (3,4)])
        >>> network = PhyloNetwork(graph)
        >>> network.nodes()
        ... [0, 1, 2, 3, 4]
        >>> network.roots()
        ... [0, 2]
        >>> network.is_phylogenetic_network()
        ... True
    
        """
        roots = list(filter(self.is_root, self.nodes()))
        roots.sort()
        return roots

    @memoize_method
    def labelled_nodes(self):
        """
        Returns the set of labelled nodes.
        
        EXAMPLE::
        
           >>> network = PhyloNetwork(eNewick="((A,B,C),1);")
           >>> network.nodes()
           ... ['_6', '_5', '_4', '_3', '_2', '_1']
           >>> network.labelled_nodes()
           ... ['_6', '_5', '_4', '_3']
           >>> list(map(network.label, network.labelled_nodes()))
           ... ['1', 'C', 'B', 'A']
           >>> network.unlabelled_nodes()
           ... ['_2', '_1']
           >>> list(map(network.label, network.unlabelled_nodes()))
           ... [None, None]
           
        """
        return list(self._labels.keys())

    @memoize_method
    def unlabelled_nodes(self):
        """
        Returns the set of unlabelled nodes.
        
        EXAMPLE::
        
           >>> network = PhyloNetwork(eNewick="((A,B,C),1);")
           >>> network.nodes()
           ... ['_6', '_5', '_4', '_3', '_2', '_1']
           >>> network.labelled_nodes()
           ... ['_6', '_5', '_4', '_3']
           >>> map(network.label, network.labelled_nodes())
           ... ['1', 'C', 'B', 'A']
           >>> network.unlabelled_nodes()
           ... ['_2', '_1']
           >>> map(network.label, network.unlabelled_nodes())
           ... [None, None]
           
        """
        return list(set(self.nodes()) - set(self.labelled_nodes()))

    @memoize_method
    def interior_nodes(self):
        """
        Returns the set of non-leaf nodes.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2),(3,4));")
            >>> network.nodes()
            ... ['_7', '_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.interior_nodes()
            ... ['_5', '_2', '_1']
            >>> map(network.is_leaf, network.interior_nodes())
            ... [False, False, False]
            
        """
        return list(set(self.nodes()) - set(self.leaves()))

    @memoize_method
    def elementary_nodes(self):
        """
        Return the set of elementary nodes.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1)E,2);")
            >>> network.nodes()
            ... ['_4', '_3', '_2', '_1']
            >>> map(network.is_elementary_node, network.nodes())
            ... [False, False, True, False]
            >>> network.label('_2')
            ... 'E'
            >>> network.elementary_nodes()
            ... '_2'
            
        """
        return list(filter(self.is_elementary_node, self.nodes()))

    @memoize_method
    def depth(self, u):
        """
        Returns the depth of u. If the node u is not from the 
        phylogenetic network, then returns None.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((((LEAF#1))),#1);")
            >>> network.nodes()
            ... ['#1', '_4', '_3', '_2', '_1']
            >>> map(network.depth, network.nodes())
            ... [1, 3, 2, 1, 0]
            >>> network.depth('non-existing node')
            ... None
            
        """
        if not u in self:
            return None
        return min([dijkstra_path_length(self, root, u) for root in self.roots()])

    @memoize_method
    def height(self, u):
        """
        Returns the height of u. If the node u is not from the 
        phylogenetic network, then returns None.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((((LEAF#1))),#1);")
            >>> network.nodes()
            ... ['#1', '_4', '_3', '_2', '_1']
            >>> map(network.height, network.nodes())
            ... [0, 1, 2, 3, 4]
            >>> network.height('non-existing node')
            ... None
            
        """
        if not u in self:
            return None
        if self.is_leaf(u):
            return 0
        else:
            return max(list(map(self.height, self.successors(u)))) + 1

    @memoize_method
    def mu(self, u):
        """
        Returns a tuple containing the number of paths from u to all 
        labelled nodes of the phylogenetic network.
        Returns None if the node u is not in the phylogenetic network.

        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((((LEAF1, LEAF2#1)), #1)INT,#1);")
            >>> network.taxa()
            ... ['INT', 'LEAF1', LEAF2]
            >>> network.roots()
            ... '_1'
            >>> network.mu('_1')
            ... array([1, 1, 3])
            >>> network.successors('_1')
            ... ['#1', '_2']
            >>> network.mu('#1') # This is LEAF2
            ... array([0, 0, 1])
            >>> network.mu('_2') # We lost the path root -> LEAF2
            ... array([1, 1, 2])
            
        """
        if u not in self:
            return None
        if self.is_leaf(u):
            mu = numpy.zeros(len(self.taxa()), int)
        else:
            mu = sum(map(self.mu, self.successors(u)))
        if self.is_labelled(u):
            pos = self.taxa().index(self.label(u))
            mu[pos] += 1
        return mu

    @memoize_method
    def mu_string(self):
        """
        Returns a string representing the mu value of all nodes, with the same order as sorted_nodes()
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((((LEAF1, LEAF2#1)), #1)INT,#1);")
            >>> network.taxa()
            ... ['INT', 'LEAF1', LEAF2]
            >>> network.sorted_nodes()
            ... ['#1', '_5', '_4', '_3', '_2', '_1']
            >>> network.mu_string()
            ... '[0 0 1]-[0 1 0]-[0 1 1]-[0 1 1]-[1 1 2]-[1 1 3]'
            >>> network.mu('#1')
            ... array([0, 0, 1])
            >>> network.mu('_1')
            ... array([1, 1, 3])
            
        """
        return '-'.join([str(self.mu(u)) for u in self.sorted_nodes()])

    @memoize_method
    def sorted_nodes(self):
        """
        Returns the set of nodes sorted with the total order over their mu value.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((((LEAF1, LEAF2#1)), #1)INT,#1);")
            >>> network.taxa()
            ... ['INT', 'LEAF1', LEAF2]
            >>> network.sorted_nodes()
            ... ['#1', '_5', '_4', '_3', '_2', '_1']
            >>> network.mu_string()
            ... '[0 0 1]-[0 1 0]-[0 1 1]-[0 1 1]-[1 1 2]-[1 1 3]'
            >>> network.nodes()
            ... ['_5', '_4', '_3', '_2', '_1', '#1']
            
        """
        nodes = self.nodes()[:]
        nodes.sort(key=cmp_to_key(lambda u, v: total_cmp(self.mu(u), self.mu(v))))
        return nodes

    def _generate_new_id(self):
        """
        For private use, it generates a new identification for every node 
        when generating the phylogenetic network.
        """

        try:
            self._lastlabel += 1
        except:
            self._lastlabel = 1
        return '_%d' % (self._lastlabel)

    # getlabel = _generate_new_id # DEPRECATED

    def _walk(self, parsed, ignore_prefix=None):
        if isinstance(parsed, pyparsing.ParseResults):

            if 'tag' in parsed:
                internal_label = '#' + str(parsed['tag'])
            else:
                internal_label = self._generate_new_id()
            if 'length' in parsed:
                length = float(parsed['length'][0])
            else:
                length = None
            logger.debug(f'adding node {internal_label}')
            self.add_node(internal_label)
            if 'label' in parsed:
                logger.debug(f'setting label {parsed["label"]} to node {internal_label}')
                self._labels[internal_label] = parsed['label']
            for child in parsed:
                walk_child = self._walk(child, ignore_prefix=ignore_prefix)
                if walk_child:
                    child_label, child_length = walk_child
                    logger.debug(f'adding arc from {internal_label} to {child_label}')
                    self.add_edge(internal_label, child_label)
                    if child_length is not None:
                        logger.debug(f'setting arc length from {internal_label} to {child_label} with length {child_length}')
                        self.edge[internal_label][child_label]['length']=child_length
            return internal_label, length

    def _from_eNewick(self, string, ignore_prefix=None):
        try:
            parsed = eNewickParser(string)[0]
        except pyparsing.ParseException:
            raise MalformedNewickException("Malformed eNewick string")
        self._walk(parsed, ignore_prefix=ignore_prefix)
        self.cache = {}

    def _eNewick_node(self, u, visited):
        if self.is_leaf(u):
            # return self._labels[u]
            return self._labels.get(u, '')
        if u in visited:
            return u
        visited.append(u)

        def length_suffix(x):
            e = self.edge[u][x]
            return ":{}".format(e['length']) if 'length' in e else ""

        children = [self._eNewick_node(x, visited) + length_suffix(x) for x in self.successors(u)]
        children.sort()
        internal = ','.join(children)
        mylabel = self.label(u) or ''
        if self.is_hybrid_node(u):
            mylabel += u
        return '(' + internal + ')' + mylabel

    def __str__(self):
        return self.eNewick()

    def __repr__(self):
        return "Phylogenetic Network with taxa [" + ",".join(map(str, self.taxa())) + "]."

    def eNewick(self):
        """
        Returns the eNewick representation of the network.
        """

        visited = []
        string = ''
        for root in self.roots():
            string += self._eNewick_node(root, visited) + ';'
        return string

    @memoize_method
    def descendant_nodes(self, u):
        """
        Returns a set with all the descendents of u.

        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((,(3,4)#1)2,#1)1;")
            >>> network.nodes()
            ... ['_5', '_4', '_3', '_2', '_1', '#1']
            >>> network.node_by_taxa('2')
            ... '_2'
            >>> network.descendant_nodes('_2')
            ... ['_5', '_4', '#1' '_3', '_2']
            >>> network.strict_descendant_nodes('_2')
            ... ['_3', '_2']
            >>> network.descendant_taxa('_2')
            ... ['4', '3', '2']
            
        """
        return sum(list(dfs_successors(self, u).values()), [])+ [u]
        # return dfs_successors(self,u)

    @memoize_method
    def descendant_taxa(self, u):
        """
        Returns a set with all the labelled descendents of u.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((,(3,4)#1)2,#1)1;")
            >>> network.nodes()
            ... ['_5', '_4', '_3', '_2', '_1', '#1']
            >>> network.node_by_taxa('2')
            ... '_2'
            >>> network.descendant_taxa('_2')
            ... ['4', '3', '2']
            >>> network.strict_descendant_taxa('_2')
            ... ['2']
            >>> network.descendant_nodes('_2')
            ... ['_5', '_4', '#1', '_3', '_2']
            
        """
        return [self.label(desc) for desc in self.descendant_nodes(u) if self.is_labelled(desc)]

    @memoize_method
    def strict_descendant_nodes(self, u):
        """
        Returns a set with all the strict descendents of u.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((,(3,4)#1)2,#1)1;")
            >>> network.nodes()
            ... ['_5', '_4', '_3', '_2', '_1', '#1']
            >>> network.node_by_taxa('2')
            ... '_2'
            >>> network.descendant_nodes('_2')
            ... ['_5', '_4', '#1', '_3', '_2']
            >>> network.strict_descendant_nodes('_2')
            ... ['_3', '_2']
            
        """
        if self.is_root(u):
            return self.descendant_nodes(u)
        pruned = copy.deepcopy(self)
        pruned.cache = {}
        pruned.remove_node(u)
        desc_pruned = []
        for root in self.roots():
            desc_pruned.extend(pruned.descendant_nodes(root))  # dfs_successors(pruned, root)
        return [desc for desc in self.descendant_nodes(u) if not desc in desc_pruned]

    @memoize_method
    def strict_descendant_taxa(self, u):
        """
        Returns a set with all the strict labelled descendents of u.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((,(3,4)#1)2,#1)1;")
            >>> network.nodes()
            ... ['_5', '_4', '_3', '_2', '_1', '#1']
            >>> network.node_by_taxa('2')
            ... '_2'
            >>> network.descendant_taxa('_2')
            ... ['4', '3', '2']
            >>> network.strict_descendant_taxa('_2')
            ... ['2']
            
        """
        return [self.label(desc) for desc in self.strict_descendant_nodes(u) if self.is_labelled(desc)]

    @memoize_method
    def ancestors(self, taxon):
        """
        Returns a set with all nodes that have a fixed descendant taxa.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((3,4)1,(5,6,7)2);")
            >>> network.ancestors('3')
            ... ['_3', '_2', '_1']
            >>> '3' in network.descendant_taxa('_2')
            ... True
            
        """

        return [u for u in self.sorted_nodes() if taxon in self.descendant_taxa(u)]

    @memoize_method
    def strict_ancestors(self, taxon):
        """
        Returns a set with all nodes that have a fixed strict descendant taxa.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick=((,(3,4)#1)2,#1)1;)
            >>> network.node_by_taxa('2')
            ... '_2'
            >>> '_2' in network.ancestors('3')
            ... True
            >>> '_2' in network.strict_ancestors('3')
            ... False
            
        """
        return [u for u in self.sorted_nodes() if taxon in self.strict_descendant_taxa(u)]

    @memoize_method
    def CSA(self, tax1, tax2):
        """
        Returns a set with the common strict ancestors of taxa1 and taxa2.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick=((,(3,4)#1)2,#1)1;)
            >>> network.CSA('3', '4')
            ... ['#1', '_1']
            >>> network.LCSA('3', '4')
            ... '#1'
            
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="(((1)), 2);")
            >>> network.CSA('1', '2')
            ... ['_2', '_1']
            >>> network.LCSA('1', '2')
            ... '_2'
            
        """

        return [u for u in self.ancestors(tax1) if
                (u in self.ancestors(tax2)) and
                ((u in self.strict_ancestors(tax1)) or
                 (u in self.strict_ancestors(tax2)))]

    @memoize_method
    def LCSA(self, tax1, tax2):
        """
        Returns a minimum of CSA(taxa1, taxa2) respect the height of nodes.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick=((,(3,4)#1)2,#1)1;)
            >>> network.CSA('3', '4')
            ... ['#1', '_1']
            >>> network.LCSA('3', '4')
            ... '#1'
            
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="(((1)), 2);")
            >>> network.CSA('1', '2')
            ... ['_2', '_1']
            >>> network.LCSA('1', '2')
            ... '_2'
            
        """
        csa = self.CSA(tax1, tax2)
        # print(self,tax1,tax2,csa)
        csa.sort(key=self.height)
        return csa[0]

    @memoize_method
    def nodal_matrix(self):
        """
        Returns a matrix containing the nodal 'distance' between all labelled nodes.
        
        EXAMPLES::
        
            >>> network = PhyloNetwork(eNewick="(((1,2), 3), 4);")
            >>> network.nodal_matrix()
            ... array([[0, 1, 2, 3],
            ...       [1, 0, 2, 3],
            ...       [1, 1, 0, 2], 
            ...       [1, 1, 1, 0])
            
        """
        n = len(self.taxa())
        matrix = numpy.zeros((n, n), int)
        dicdist = all_pairs_shortest_path_length(self)
        for i in range(n):
            ti = self.taxa()[i]
            for j in range(i, n):
                tj = self.taxa()[j]
                lcsa = self.LCSA(ti, tj)
                matrix[i, j] = dicdist[lcsa][self.node_by_taxa(ti)]
                matrix[j, i] = dicdist[lcsa][self.node_by_taxa(tj)]
        return matrix

    def nodal_area(self):
        """
        Returns the sum of all elements of the nodal matrix.
        
        EXAMPLES::
        
            >>> network = PhyloNetwork(eNewick="(((1,2), 3), 4);")
            >>> network.nodal_matrix()
            ... array([[0, 1, 2, 3],
            ...       [1, 0, 2, 3],
            ...       [1, 1, 0, 2], 
            ...       [1, 1, 1, 0])
            >>> network.nodal_area()
            ... 19
            
        """
        mat = self.nodal_matrix()
        # mat=mat+mat.transpose()
        return sum(abs(mat.flatten()))

    @memoize_method
    def cophenetic_matrix(self):
        """
        Returns a matrix with the cophenetic coeficient of labelled nodes.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="(((1,2), 3), 4);")
            >>> network.cophenetic_matrix()
            ... array([[3, 2, 1, 0],
            ...       [0, 3, 1, 0],
            ...       [0, 0, 2, 0], 
            ...       [0, 0, 0, 1])
            
        """
        n = len(self.taxa())
        matrix = numpy.zeros((n, n), int)
        for i in range(n):
            ti = self.taxa()[i]
            for j in range(i, n):
                tj = self.taxa()[j]
                lcsa = self.LCSA(ti, tj)
                matrix[i, j] = self.depth(lcsa)
        return matrix

    def common_taxa(self, net2):
        """
        Returns common taxa between self and net2.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2), 3)4;")
            >>> network2 = PhyloNewtwork(eNewick="(1,(4,5)2);")
            >>> network.common_taxa(network2)
            ... ['1', '2', '4']
            >>> network.common_taxa_leaves(network2)
            ... ['1']
            
        """

        common = []
        taxa1 = self.taxa()
        taxa2 = net2.taxa()
        for taxon in taxa1:
            if taxon in taxa2:
                common.append(taxon)
        return common

    def common_taxa_leaves(self, net2):
        """
        Returns common taxa between self and net2 that are leaves on both networks.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2), 3)4;")
            >>> network2 = PhyloNewtwork(eNewick="(1,(4,5)2);")
            >>> network.common_taxa(network2)
            ... ['1', '2', '4']
            >>> network.common_taxa_leaves(network2)
            ... ['1']
            
        """

        common = []
        taxa1 = list([l for l in self.taxa() if self.is_leaf(self.node_by_taxa(l))])
        taxa2 = list([l for l in net2.taxa() if net2.is_leaf(net2.node_by_taxa(l))])
        for taxon in taxa1:
            if taxon in taxa2:
                common.append(taxon)
        return common

    def topological_restriction(self, subtaxa, nested=True):
        """
        Returns a minimal subnetwork of self such that it containts a fixed subtaxa.
        If nested=False then only taxa on leaves will be considered.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2), (3,4));")
            >>> network.taxa()
            ... ['1', '2', '3', '4']
            >>> network.roots()
            ... ['_1']
            >>> subnetwork = network.topological_restriction(['1', '2'])
            >>> subnetwork.taxa()
            ... ['1', '2']
            >>> subnetwork.nodes() # '_1' should not be
            ... ['_4', '_3', '_2']
            >>> subnetwork.roots()
            ... ['_2']
            >>> subnetwork = network.topological_restriction(['1', '3'])
            >>> '_1' in subnetwork.roots()
            ... True
            
        """
        restricted = copy.deepcopy(self)
        for node in restricted.labelled_nodes():
            if restricted.label(node) not in subtaxa:
                # Delete taxa not in subtaxa
                del restricted._labels[node]
            elif not nested and not restricted.is_leaf(node):
                # if not nested, all taxa (even if they are in subtaxa) 
                # on internal nodes should be removed.
                del restricted._labels[node]
        restricted.cache = {}

        candidates_to_delete = restricted.leaves()
        while len(candidates_to_delete) > 0:
            node = candidates_to_delete.pop()
            if restricted.is_leaf(node) and not restricted.is_labelled(node):
                candidates_to_delete[0:0] = restricted.predecessors(node)
                restricted.remove_node(node)
                restricted.cache = {}

        for u in restricted.nodes():
            if restricted.is_elementary_node(u) and \
                    not restricted.is_labelled(u):
                for parent in restricted.predecessors(u):
                    for child in restricted.successors(u):
                        restricted.add_edge(parent, child)
                restricted.remove_node(u)
        restricted.cache = {}
        return restricted

    @memoize_method
    def has_nested_taxa(self):
        """
        Returns True is an internal node is labelled. False otherwise.
        """

        for node in self.labelled_nodes():
            if not self.is_leaf(node):
                return True
        return False

    def matching_representation(self):
        try:
            return self._matching_representation
        except:
            self._matching_representation = {}
            for u in self.nodes():
                self._matching_representation[u] = 0
            for u in self.leaves():
                pos = self.taxa().index(self.label(u))
                self._matching_representation[u] = pos + 1
            h = 1
            i = len(self.leaves()) + 1
            thislevel = [u for u in self.nodes() if self.height(u) == h]
            while thislevel:
                minims = {}
                for u in thislevel:
                    minims[u] = min([self._matching_representation[v] for \
                                     v in self.successors(u)])
                thislevel.sort(key = cmp_to_key(lambda x, y: minims[x] - minims[y]))
                for k in range(len(thislevel)):
                    self._matching_representation[thislevel[k]] = i
                    i += 1
                h += 1
                thislevel = [u for u in self.nodes() if self.height(u) == h]
            return self._matching_representation

    def matching_permutation(self):
        try:
            return self._matching_permutation
        except:
            permutation = {}
            representation = self.matching_representation()
            for u in self.nodes():
                children = self.successors(u)
                children.sort(key=cmp_to_key(lambda u, v: representation[u] - representation[v]))
                for i in range(len(children)):
                    permutation[representation[children[i]]] = representation[children[(i + 1) % len(children)]]
            self._matching_permutation = permutations.Permutation(permutation)
            return self._matching_permutation

    def cluster(self, u):
        """
        Returns the cluster of u in the network, None if the node is not in the network.
        """

        if not u in self:
            return None

        cl = []
        dictio = single_source_shortest_path_length(self, u)
        for node in list(dictio.keys()):
            if self.label(node):
                cl.append(self.label(node))
        cl.sort()
        cl = tuple(cl)
        return cl

    def cluster_representation(self):
        """
        Returns the cluster of all nodes in the network.
        """

        cls = list(map(self.cluster, self.nodes()))
        cls.sort()
        return cls

    def nested_label(self, node):
        """
        Returns a string representation of descendants of u. Very useful to identify where is a node located in the network.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((2,3), (4,(5,6)))1;")
            >>> network.nodes()
            ... ['_9', '_8', '_7', '_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.nested_label('_7') # Where is node '_7' in the network?
            ... '{5,6}' # '_7' is the parent of leaves '5' and '6'.
            >>> network.nested_label('_6')
            ... '4' # '_6' is the leaf '4'
            >>> network.nested_label('_5')
            ... '{4,{5,6}}' # '_5' is the parent of leaf '4' and node '_7'
        """

        try:
            return self._nested_label[node]
        except:
            if not hasattr(self, '_nested_label'):
                self._nested_label = {}
            if self.is_leaf(node):
                return self.label(node)
            else:
                desc_labels = [self.nested_label(u) for u in self.successors(node)]
                desc_labels.sort()
                if desc_labels is None:
                    # If a leaf doesn't have a label
                    self._nested_label[node] = "{}"
                    return self.nested_label[node]
                self._nested_label[node] = '{' + ','.join(desc_labels) + '}'
                return self._nested_label[node]

    def nested_label_representation(self, multiset=False):
        """
        Returns the nested label of all nodes in the network. If the optional parameter
        is True, then the result is returned as a Counter object; otherwise as a set
        object.
        
        EXAMPLE::
        
            >>> network = PhyloNetwork(eNewick="((1,2), (3,4));")
            >>> network.nodes()
            ... ['_7', '_6', '_5', '_4', '_3', '_2', '_1']
            >>> network.nested_label_representation()
            ... set(['{3,4}', '1', '3', '2', '4', '{1,2}', '{{1,2},{3,4}}'])
            
        """
        nls = [self.nested_label(u) for u in self.nodes()]
        if multiset:
            return Counter(nls)
        else:
            return set(nls)

    def draw(self):
        try:
            pydot_layout = nx.nx_pydot.pydot_layout
            draw_networkx = nx.draw_networkx()
            pos = nx.pydot_layout(self, prog='dot')
            nx.draw_networkx(self, pos, labels=self._labels, edgelist=[])
            nx.draw_networkx_edges(self, pos)
        except:
            pass


def _get_chunck(string):
    if string[0] in '(),;':
        return string[0], string[1:]
    i = 0
    while string[i] not in '(),;':
        i += 1
    return string[0:i], string[i:]


class PhyloTree(PhyloNetwork):
    def _from_eNewick(self, string, ignore_prefix=None):
        substr = string
        current = self._generate_new_id()
        self.add_node(current)
        finished = False
        while not finished:
            chunk, substr = _get_chunck(substr)
            # print(chunk,substr)
            if chunk == '(':
                child = self._generate_new_id()
                self.add_edge(current, child)
                current = child
                continue
            if chunk == ')':
                current = self.predecessors(current)[0]
                continue
            if chunk == ',':
                parent = self.predecessors(current)[0]
                child = self._generate_new_id()
                self.add_edge(parent, child)
                current = child
            if chunk == ';':
                finished = True
                continue
            else:
                # chunk is label:dist
                if ':' in chunk:
                    spl = chunk.split(':')
                    label = spl[0]
                    length = spl[1]
                else:
                    label = chunk
                    length = None
                if label:
                    self._labels[current] = label
                if length:
                    try:
                        length_value = int(length)
                    except:
                        length_value = float(length)
                    parent = self.predecessors(current)[0]
                    self.edge[parent][current]['length'] = length_value



                #    @memoize_method
                #    def cluster(self,u):
                #        if self.is_leaf(u):
                #            return set([self.label(u)])
                #        cl = set()
                #        for v in self.successors(u):
                #            cl = cl | self.cluster(v)
                #        return cl
                #
                #    @memoize_method
                #    def depth(self,u):
                #        if self.is_root(u):
                #            return 0
                #        return 1+self.depth(self.predecessors(u)[0])

                #    @memoize_method
                #    def nodal_matrix(self):
                #        taxa = self.taxa()
                #        n = len(taxa)
                #        taxa_map = {taxa[i]:i for i in range(n)}
                #        matrix=numpy.zeros((n,n),int)
                #        for u in self.interior_nodes():
                #            descendant_clusters = [self.cluster(v) for v in self.successors(u)]
                #            pairs = combinations(descendant_clusters,2)
                #            for pair in pairs:
                #                cl1 = pair[0]
                #                cl2 = pair[1]
                #                for taxa_i in cl1:
                #                    for taxa_j in cl2:
                #                        i = taxa_map[taxa_i]
                #                        j = taxa_map[taxa_j]
                #                        li = self.depth(self.node_by_taxa(taxa_i)) - self.depth(u)
                #                        lj = self.depth(self.node_by_taxa(taxa_j)) - self.depth(u)
                #                        matrix[i,j] = li
                #                        matrix[j,i] = lj
                #        return matrix

    @memoize_method
    def cophenetic_matrix(self):
        taxa = self.taxa()
        n = len(taxa)
        taxa_map = {taxa[i]: i for i in range(n)}
        matrix = numpy.zeros((n, n), int)
        for u in self.interior_nodes():
            descendant_clusters = [self.cluster(v) for v in self.successors(u)]
            pairs = combinations(descendant_clusters, 2)
            for pair in pairs:
                cl1 = pair[0]
                cl2 = pair[1]
                for taxa_i in cl1:
                    for taxa_j in cl2:
                        i = taxa_map[taxa_i]
                        j = taxa_map[taxa_j]
                        matrix[min(i, j), max(i, j)] = self.depth(u)
        for u in self.leaves():
            i = taxa_map[self.label(u)]
            matrix[i, i] = self.depth(u)
        return matrix

    def remove_elementary_nodes(self):
        elementary_nodes = [u for u in self.nodes() if (self.out_degree(u) == 1 and not self.is_labelled(u))]
        for u in elementary_nodes:
            child = self.successors(u)[0]
            for v in self.predecessors(u):
                self.add_edge(v, child)
            self.remove_node(u)

    @memoize_method
    def cluster(self, u):
        """
        Returns the cluster of u in the tree, None if the node is not in the network.
        """

        if not u in self:
            return None

        if self.is_labelled(u):
            return set([self.label(u)])
        else:
            tmp = set()
            for child in self.successors(u):
                tmp = tmp.union(self.cluster(child))
            return tmp

    @memoize_method
    def depth(self, u):
        """
        Returns the depth of u. If node u is not from the 
        tree, then returns None.
        """
        if not u in self:
            return None

        if self.is_root(u):
            return 0
        else:
            return 1 + self.depth(self.predecessors(u)[0])

    @memoize_method
    def _compute_LCA(self, i, j):
        if i == j:
            return i
        if self.is_root(i):
            return i
        parent = self.predecessors(i)[0]
        if j in self.descendant_nodes(parent):
            return parent
        return self._compute_LCA(parent, j)

    @memoize_method
    def LCA(self, tax1, tax2):
        """
        Returns the minimum of CSA(taxa1, taxa2) respect the height of nodes.
            
        EXAMPLE::
        
            >>> network = PhyloTree(eNewick="(((1)), 2);")
            >>> network.CSA('1', '2')
            ... ['_2', '_1']
            >>> network.LCSA('1', '2')
            ... '_2'
            
        """

        node1 = self.node_by_taxa(tax1)
        node2 = self.node_by_taxa(tax2)
        return self._compute_LCA(node1, node2)

    LCSA = LCA

    @memoize_method
    def nodal_matrix(self):
        """
        Returns a matrix containing the nodal 'distance' between all labelled nodes.
        
        EXAMPLES::
        
            >>> network = PhyloTree(eNewick="(((1,2), 3), 4);")
            >>> network.nodal_matrix()
            ... array([[0, 1, 2, 3],
            ...       [1, 0, 2, 3],
            ...       [1, 1, 0, 2], 
            ...       [1, 1, 1, 0])
            
        """
        n = len(self.taxa())
        matrix = numpy.zeros((n, n), int)
        for i in range(n):
            ti = self.taxa()[i]
            for j in range(i, n):
                tj = self.taxa()[j]
                lca = self.LCA(ti, tj)
                matrix[i, j] = self.depth(self.node_by_taxa(ti)) - self.depth(lca)
                matrix[j, i] = self.depth(self.node_by_taxa(tj)) - self.depth(lca)
        return matrix


class LGTPhyloNetwork(PhyloNetwork):
    def _walk(self, parsed, ignore_prefix=None):
        if isinstance(parsed, pyparsing.ParseResults):
            secondary = False
            if 'tag' in parsed:
                internal_label = '#' + str(parsed['tag'])
                if 'type' in parsed and parsed['type'] == 'LGT':
                    secondary = True
            else:
                internal_label = self._generate_new_id()
            if 'length' in parsed:
                pass
            self.add_node(internal_label)
            if 'label' in parsed:
                self._labels[internal_label] = parsed['label']
            for child in parsed:
                inner_walk = self._walk(child, ignore_prefix=ignore_prefix)
                if inner_walk:
                    child_label = inner_walk[0]
                    secondary_edge = inner_walk[1]
                    self.add_edge(internal_label, child_label, attr_dict={'secondary': secondary_edge})
            return [internal_label, secondary]

    def principal_edges(self):
        return [(edge[0], edge[1]) for edge in self.edges(data=True) if
                not 'secondary' in edge[2] or edge[2]['secondary'] == False]

    def secondary_edges(self):
        return [(edge[0], edge[1]) for edge in self.edges(data=True) if
                'secondary' in edge[2] and edge[2]['secondary'] == True]

    def principal_subtree(self, simplify=True):
        tree = PhyloTree()
        tree.add_edges_from(self.principal_edges())
        tree._labels = self._labels
        if simplify:
            tree.remove_elementary_nodes()
        return tree

    def secondary_subtrees(self, simplify=True):
        subtrees = []
        for edge in self.secondary_edges():
            tree = PhyloTree()
            edges = [e for e in self.principal_edges() if e[1] != edge[1]]
            edges.append(edge)
            tree.add_edges_from(edges)
            tree._labels = self._labels
            if simplify:
                tree.remove_elementary_nodes()
            subtrees.append(tree)
        return subtrees

    def draw(self):
        try:
            pydot_layout = nx.nx_pydot.pydot_layout
            draw_networkx = nx.draw_networkx()
            pos = nx.pydot_layout(self, prog='dot')
            nx.draw_networkx(self, pos, labels=self._labels, edgelist=[])
            nx.draw_networkx_edges(self, pos, edgelist=self.principal_edges())
            nx.draw_networkx_edges(self, pos, edgelist=self.secondary_edges(), style='dashed')
        except:
            pass
