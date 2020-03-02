from collections import Counter
from itertools import combinations_with_replacement
from networkx import union
from functools import lru_cache
from .utils import product as prod
from math import factorial
from . import PhylogeneticTree, NetworkShape

try:
    from functools import cached_property
except ImportError:
    from cached_property import cached_property
import logging
logger = logging.getLogger(__name__)


class TreeShape(NetworkShape):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def enroot(cls, t1, t2):
        r1 = f'1{t1.root}'
        r2 = f'2{t2.root}'
        result = union(t1, t2, rename=('1', '2'))
        result.add_edge('0', r1)
        result.add_edge('0', r2)
        result.__class__ = TreeShape
        return result

    @classmethod
    def leaf(cls):
        result = TreeShape()
        result.add_node('0')
        return result

    def newick(self, u):
        if self.is_leaf(u):
            return '*'
        return '(' + ','.join([self.newick(v) for v in self.successors(u)]) + ')'

    def __str__(self):
        return self.newick(self.root) + ';'

    @lru_cache()
    def canonical_label(self, u):
        if self.is_leaf(u):
            return '10'
        return '1' + ''.join(sorted([self.canonical_label(v) for v in self.successors(u)])) + '0'

    def number_of_symmetries_at_node(self, u):
        if self.is_leaf(u):
            return 1
        labels = [self.canonical_label(v) for v in self.successors(u)]
        c = Counter(labels)
        return prod(map(factorial, c.values()))

    def number_of_symmetries(self):
        return prod([self.number_of_symmetries_at_node(u) for u in self.nodes])

    def number_of_phylogenetic_trees(self):
        n = len(self.leaves)
        return factorial(n) // self.number_of_symmetries()

    def sample_phylogenetic_tree(self):
        net = self.copy()
        net.__class__ = PhylogeneticTree
        for (i, leaf) in enumerate(net.leaves):
            net.nodes[leaf]['label'] = str(i+1)
        return net


@lru_cache()
def maximal_balanced_tree_shape(n):
    if n == 1:
        return TreeShape.leaf()
    if n % 2 == 0:
        n1 = n // 2
        n2 = n // 2
    else:
        n1 = n // 2
        n2 = n1 + 1
    t1 = maximal_balanced_tree_shape(n1)
    t2 = maximal_balanced_tree_shape(n2)
    return TreeShape.enroot(t1, t2)


@lru_cache()
def caterpillar_tree_shape(n):
    if n == 1:
        return TreeShape.leaf()
    return TreeShape.enroot(TreeShape.leaf(), caterpillar_tree_shape(n-1))


@lru_cache()
def all_binary_tree_shapes(n):
    if n == 1:
        return [TreeShape.leaf()]
    if n % 2 == 0:
        maxi = n//2 - 1
    else:
        maxi = n//2
    results = []
    for i in range(1, maxi + 1):
        j = n-i
        for t1 in all_binary_tree_shapes(i):
            for t2 in all_binary_tree_shapes(j):
                results.append(TreeShape.enroot(t1, t2))
    if n % 2 == 0:
        i = n//2
        for (t1, t2) in combinations_with_replacement(all_binary_tree_shapes(i), 2):
            results.append(TreeShape.enroot(t1, t2))
    return results
