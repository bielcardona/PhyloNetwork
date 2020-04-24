"""
Main module of the package. Defines classes and methods for working with phylogenetic trees and networks.
"""

import networkx as nx
import numpy
from networkx import DiGraph
import pyparsing
from .eNewick import eNewickParser
from .exceptions import *
from itertools import combinations

from cached_property import cached_property

import logging
logger = logging.getLogger(__name__)


class NetworkShape(DiGraph):
    """
    Class for unlabeled phylogenetic networks/trees.

    In this class we only consider properties that depend on the topology of the network, without a labeling.

    We implicitly assume:

    * There is a single root
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._last_id = 0
        self._last_reticulation_id = 0

    def _generate_new_id(self):
        """
        For private use, generates new id for nodes.
        """
        self._last_id += 1
        while f'_{self._last_id}' in self.nodes:
            self._last_id += 1
        return f'_{self._last_id}'

    def _generate_new_reticulation_id(self):
        """
        For private use, generates new id for nodes.
        """
        self._last_reticulation_id += 1
        while f'#{self._last_reticulation_id}' in self.nodes:
            self._last_reticulation_id += 1
        return f'#{self._last_reticulation_id}'

    def new_node(self):
        """Creates a node as a string of the form ``_n`` with unique n, and returns it"""
        u = self._generate_new_id()
        self.add_node(u)
        return u

    def new_reticulation_node(self):
        """Creates a node as a string of the form ``#n`` with unique n, and returns it"""
        u = self._generate_new_reticulation_id()
        self.add_node(u)
        return u

    def is_leaf(self, u):
        """Tests if `u` is a leaf of the network"""
        return self.out_degree(u) == 0  # noqa

    @cached_property
    def leaves(self):
        """The set of all leaves of the network"""
        return set([u for u in self.nodes if self.is_leaf(u)])

    @cached_property
    def interior_nodes(self):
        """The set of all interior (not leaves) nodes of the network"""
        return set([u for u in self.nodes if not self.is_leaf(u)])

    def is_reticulation(self, u):
        """Tests if node `u` is a reticulation (indeg >= 2)"""
        return self.in_degree(u) > 1  # noqa

    @cached_property
    def reticulations(self):
        """The set of all reticulations"""
        return set([u for u in self.nodes if self.is_reticulation(u)])

    def is_tree_node(self, u):
        """Tests if node `u` is a tree node (indeg <= 1)"""
        return self.in_degree(u) <= 1  # noqa

    @cached_property
    def tree_nodes(self):
        """The set of all tree nodes"""
        return set([u for u in self.nodes if self.is_tree_node(u)])

    def is_root(self, u):
        """Tests if `u` is a root"""
        return self.in_degree(u) == 0  # noqa

    @cached_property
    def roots(self):
        """The set of all roots of the network"""
        return set([u for u in self.nodes if self.is_root(u)])

    @cached_property
    def root(self):
        """Root of the network. Raises exception if no root or >1 roots"""
        if len(self.roots) != 1:
            raise NotAPhylogeneticNetwork
        for u in self.roots:
            return u

    def is_elementary(self, u):
        """Tests if `u` is elementary"""
        return self.in_degree(u) <= 1 and self.out_degree(u) == 1  # noqa

    @cached_property
    def elementary_nodes(self):
        """The set of elementary nodes"""
        return set([u for u in self.nodes if self.is_elementary(u)])

    @cached_property
    def bottom_to_top_nodes(self):
        """List of nodes ordered from bottom to top (last node is the root)"""
        return list(reversed(list(nx.topological_sort(self))))

    @cached_property
    def top_to_bottom_nodes(self):
        """List of nodes ordered from top to bottom (first node is the root)"""
        return list(nx.topological_sort(self))

    @cached_property
    def depths(self):
        """Dict associating each node to its depth (distance to root)"""
        return nx.single_source_shortest_path_length(self, self.root)

    @cached_property
    def heights(self):
        """Dict associating each node to its height (max. distance to a leaf)"""
        logger.info('computing heights')
        _heights = {}
        for u in self.bottom_to_top_nodes:
            if self.is_leaf(u):
                _heights[u] = 0
            else:
                _heights[u] = 1 + max([_heights[v] for v in self.successors(u)])
        return _heights

    def clear_cache(self):
        """Clears the cache of all computed properties"""
        names = dir(self.__class__)
        for name in names:
            attr = getattr(self.__class__, name)
            if attr.__class__ is cached_property:
                logger.info(f'found chached property {name}...')
                try:
                    delattr(self, name)
                    logger.info(f'deleted!')
                except:
                    logger.info(f'could not delete!')

    def remove_node_and_reconnect(self, u, clear_cache=False):
        """Removes node `u` and connects each of its parents to each of its children"""
        for parent in self.predecessors(u):
            for child in self.successors(u):
                self.add_edge(parent, child)
        self.remove_node(u)
        if clear_cache:
            self.clear_cache()

    def remove_elementary_nodes(self):
        """Removes all elementary nodes of the network"""
        while self.elementary_nodes:
            for u in self.elementary_nodes:
                self.remove_node_and_reconnect(u)
            self.clear_cache()

    def make_elementary_node(self, u, v):
        """Adds an elementary node between `u` and `v`. Returns the added node."""
        self.remove_edge(u, v)
        w = self.new_node()
        self.add_edges_from([(u, w), (w, v)])
        return w

    def level(self):
        """Returns the level of the network"""
        net = nx.Graph(self)
        bcc = nx.biconnected_components(net)
        no_trivial_bcc = [x for x in bcc if len(x) > 2]
        if len(no_trivial_bcc) == 0:
            return 0
        return max([len(([x for x in c if self.is_reticulation(x)])) for c in no_trivial_bcc])


class PhylogeneticNetwork(NetworkShape):
    """
    Class for (labeled) phylogenetic networks.

    We assume:

    * All leaves are labeled
    * Internal nodes _may_ be labeled
    * All labels are different

    If `eNewick` is given, it will be populated.

    The (e)Newick features implemented are:

    * Reticulation events [CRV08]_
    * Additional data on nodes and arcs. Format:

    ::

        [&node_metadata]:[&arc_metadata]arc_length

    """

    def __init__(self, eNewick=None, **kwargs):  # noqa
        """
        Initializes a PhylogeneticNetwork instance. If `eNewick` is given, it will be populated.

        The (e)Newick features implemented are:

        * Reticulation events [CRV08]
        * Node and arc metadata
        * Branch lenghts
        """
        super().__init__(kwargs)
        self._last_id = 0
        if eNewick is not None:
            self._from_eNewick(eNewick)

    @staticmethod
    def _process_metadata(string):
        """For private use. Processes the metadata stored in `string`"""
        result = {}
        if string[0] != '&':
            return
        string = string[1:]
        pairs = string.split(',')
        for pair in pairs:
            try:
                k, v = pair.split('=')
                result[k] = v
            except:
                pass
        return result

    def _walk(self, parsed, ignore_prefix=None):
        """For private use. Recursively processes a parsed eNewick string"""
        if isinstance(parsed, pyparsing.ParseResults):
            if 'tag' in parsed:
                node_id = '#' + str(parsed['tag'])
            else:
                node_id = self._generate_new_id()
            length = parsed.get('length', None)
            arc_metadata = parsed.get('arc_metadata', None)
            if arc_metadata is not None:
                arc_metadata = self._process_metadata(arc_metadata[0])
            logger.debug(f'adding node {node_id}')
            self.add_node(node_id)
            if 'node_metadata' in parsed:
                node_metadata = self._process_metadata(parsed['node_metadata'][0])
                logger.debug(f"adding metadata {node_metadata} to node {node_id}")
                self.nodes[node_id]['metadata'] = node_metadata
            if 'label' in parsed:
                logger.debug(f'setting label {parsed["label"]} to node {node_id}')
                self.nodes[node_id]['label'] = parsed['label']
            for child in parsed:
                walk_child = self._walk(child, ignore_prefix=ignore_prefix)
                if walk_child:
                    child_label, child_length, child_arc_metadata = walk_child
                    logger.debug(f'adding arc from {node_id} to {child_label}')
                    self.add_edge(node_id, child_label)
                    if child_length is not None:
                        logger.debug(
                            f'setting arc length from {node_id} to {child_label} with length {child_length}')
                        self.edges[node_id, child_label]['length'] = child_length
                    if child_arc_metadata is not None:
                        logger.debug(
                            f'setting arc metadata from {node_id} to {child_label} with {child_arc_metadata}')
                        self.edges[node_id, child_label]['metadata'] = child_arc_metadata

            return node_id, length, arc_metadata

    def _from_eNewick(self, string, ignore_prefix=None):
        """For private use, builds network from a given eNewick string"""
        try:
            parsed = eNewickParser(string)[0]
        except pyparsing.ParseException:
            raise MalformedNewickException("Malformed eNewick string")
        self._walk(parsed, ignore_prefix=ignore_prefix)
        self.cache = {}

    @cached_property
    def labeling_dict(self):
        """Dict that maps nodes to labels"""
        return nx.get_node_attributes(self, 'label')

    @property
    def taxa(self):
        """Set of labels used for nodes"""
        return set(self.labeling_dict.values())

    @property
    def ordered_taxa(self):
        """Ordered list of taxa"""
        return sorted(self.taxa)

    def is_labeled(self, u):
        """Tests if `u` is labeled"""
        return 'label' in self.nodes[u]

    def label(self, u):
        """Returns the label of `u` (or None if it is not labeled)"""
        return self.labeling_dict.get(u, None)

    @cached_property
    def labeled_nodes(self):
        """Set of all labeled nodes"""
        return set([u for u in self.nodes if self.is_labeled(u)])

    @cached_property
    def unlabeled_nodes(self):
        """Set of all unlabeled nodes"""
        return set([u for u in self.nodes if not self.is_labeled(u)])

    @cached_property
    def lookup_label_dict(self):
        """Dict that maps each label to the corresponding node"""
        return {v: k for k, v in self.labeling_dict.items()}

    def node_with_label(self, taxon):
        """Returns the node label by `taxon` (or None if no such label)"""
        return self.lookup_label_dict.get(taxon, None)

    def is_unlabeled_elementary_or_leaf(self, u):
        """Tests if u is elementary (not labeled and (d_in,d_out)=(1,1),(0,1) or (1,0)"""
        return (self.is_elementary(u) or self.is_leaf(u)) and not self.is_labeled(u)

    def _length_suffix(self, u, x):
        """For internal use. Returns a string representing the length of the arc ux (or empty if not set)"""
        e = self.edges[u, x]
        return f":{e['length']}" if 'length' in e else ""

    def _eNewick_node_dict(self):
        """For internal use. Returns a dict that maps each node to its associated eNewick string"""
        nodes = self.bottom_to_top_nodes
        result = {}
        visited = set()
        for u in nodes:
            if self.is_leaf(u):
                result[u] = self.nodes[u].get('label', '')
            else:
                children = []
                for x in self.successors(u):
                    if x in visited:
                        children.append(x)
                    else:
                        visited.add(x)
                        children.append(result[x])
                children.sort()
                internal = ','.join(children)
                mylabel = self.nodes[u].get('label', '')
                if self.is_reticulation(u):
                    mylabel += u
                result[u] = f"({internal}){mylabel}"
        return result

    def eNewick(self):
        """Returns the eNewick representation of the network."""
        return self._eNewick_node_dict()[self.root] + ';'

    def __str__(self):
        return self.eNewick()

    def __repr__(self):
        return f"PhylogeneticNetwork(eNewick = {self.eNewick()})"

    @cached_property
    def mu_dict(self):
        """Dict that maps each node to its mu-vector"""
        mus = {}
        for u in self.bottom_to_top_nodes:
            if self.is_leaf(u):
                mus[u] = numpy.zeros(len(self.ordered_taxa), int)
            else:
                mus[u] = sum([mus[v] for v in self.successors(u)])
            if self.is_labeled(u):
                pos = self.ordered_taxa.index(self.label(u))
                mus[u][pos] += 1
        return mus

    def mu(self, u):
        """Returns the mu-vector of `u`"""
        return self.mu_dict[u]

    @cached_property
    def cluster_dict(self):
        """Dict that maps each node to its cluster"""
        cls = {}
        n = len(self.ordered_taxa)
        for u in self.bottom_to_top_nodes:
            if self.is_leaf(u):
                cls[u] = {self.label(u)}
            else:
                cldesc = [cls[v] for v in self.successors(u)]
                cls[u] = set.union(*cldesc)
                if self.is_labeled(u):
                    cls[u] |= {self.label(u)}
        return cls

    def cluster(self, u):
        """Returns the cluster of `u`"""
        return self.cluster_dict[u]

    @cached_property
    def nested_label_dict(self):
        """Dict that maps each node to its nested label representation"""
        result = {}
        for node in self.bottom_to_top_nodes:
            if self.is_leaf(node):
                result[node] = self.label(node)
            else:
                result[node] = '{' + ','.join(sorted([result[u] for u in self.successors(node)])) + '}'
        return result

    def nested_label(self, u):
        """Returns the nested label representation of `u`"""
        return self.nested_label_dict[u]

    def remove_taxa(self, taxa, prune_network=True):
        """Removes the given `taxa` from the network and, if `prune_network`, eliminates elementary nodes"""
        for label in taxa:
            u = self.node_with_label(label)
            if u:
                del self.nodes[u]['label']
        if prune_network:
            self.remove_elementary_nodes()
            self.clear_cache()

    def topological_restriction(self, taxa):
        """Returns the topological restriction of the network to the subset `taxa` of labels"""
        net = self.copy()
        net.remove_taxa(self.taxa - taxa)
        return net

    def has_nested_taxa(self):
        """Tests if the network has labeled internal nodes"""
        for u in self.interior_nodes:
            if self.label(u):
                return True
        return False

    def is_tree_child(self):
        """Tests if the network is tree-child"""
        for u in self.interior_nodes:
            if not any([self.is_tree_node(v) for v in self.successors(u)]):
                return False
        return True

    def draw(self):
        """Plots the network with labels"""
        # requires matplotlib and pygraphviz
        import matplotlib.pyplot as plt
        from networkx.drawing.nx_agraph import graphviz_layout
        pos = graphviz_layout(self, prog='dot')
        #nx.draw(self, pos)
        nx.draw_networkx_nodes(self, pos, self.tree_nodes, node_size=200, node_color='#57f542')
        nx.draw_networkx_nodes(self, pos, self.reticulations, node_size=150, node_shape='s', node_color='#57f542')
        nx.draw_networkx_edges(self, pos)
        nx.draw_networkx_labels(self, pos, labels=self.labeling_dict)
        plt.show()


class PhylogeneticTree(PhylogeneticNetwork):
    """
    Class for working with phylogenetic trees.

    Implements some methods that are not defined for networks (or are easier to implement on trees).

    If `Newick` is given, it will be populated.
    """

    def __init__(self, Newick=None, **kwargs):
        super().__init__(eNewick=Newick, **kwargs)

    @cached_property
    def LCA_dict(self):
        """Dict that associate to each pair of taxons its least common ancestor"""
        result = {}
        for u in self.top_to_bottom_nodes:
            if self.is_leaf(u):
                l = self.label(u)
                result[(l, l)] = u
            for children in self.successors(u):
                for (ch1, ch2) in combinations(children, 2):
                    for v1 in self.cluster(ch1):
                        for v2 in self.cluster(ch2):
                            result[(v1, v2)] = u
        return result

    def LCA(self, u, v):
        """Returns the LCA of `u` and `v`"""
        return self.LCA_dict[u, v]

    def nodal_matrix(self):
        """Not implemented"""
        raise exceptions.NotImplemented

    def cophenetic_matrix(self):
        """Not implemented"""
        raise exceptions.NotImplemented

    def __repr__(self):
        return f"PhylogeneticTree(eNewick = {self.eNewick()})"


class LGTNetwork(PhylogeneticNetwork):
    """Class for representing LGT networks (Not yet implemented)"""

    def __init__(self, eNewick=None, **kwargs):
        super().__init__(eNewick=eNewick, **kwargs)

    def _walk(self, parsed, ignore_prefix=None):
        """Not implemented"""
        raise exceptions.NotImplemented

    def principal_edges(self):
        """Not implemented"""
        raise exceptions.NotImplemented

    def secondary_edges(self):
        """Not implemented"""
        raise exceptions.NotImplemented

    def principal_subtree(self, simplify=True):
        """Not implemented"""
        raise exceptions.NotImplemented

    def secondary_subtrees(self, simplify=True):
        """Not implemented"""
        raise exceptions.NotImplemented

    def draw(self):
        """Not implemented"""
        raise exceptions.NotImplemented


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
