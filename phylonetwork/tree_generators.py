"""
Module that generates phylogenetic trees either sequentially or randomly.
The implemented random models are uniform and Yule distributions
"""

from . import PhylogeneticTree as PhyloTree
from .tree_operations import push_and_hang, hold_and_hang, push_and_label, hold_and_label
from .utils import random_weighted, product
from random import randint, shuffle, choice
from functools import lru_cache
import random


# Generation using cached trees
all_binary_trees_cached = {}

def all_binary_trees_by_num_taxa(num_taxa):
    if num_taxa in all_binary_trees_cached:
        return all_binary_trees_cached[num_taxa]
    if num_taxa == 1:
        all_binary_trees_cached[0] = [tree_with_single_taxon('0')]
        return all_binary_trees_cached[0]
    taxon = str(num_taxa - 1)
    trees = []
    parents = all_binary_trees_by_num_taxa(num_taxa - 1)
    for parent in parents:
        for u in parent.nodes:
            newtree = push_and_hang(parent, u, taxon)
            _ = newtree.lookup_label_dict
            trees.append(newtree)
    all_binary_trees_cached[num_taxa] = trees
    return trees

def relabel(tree, lookup_dict, taxa):
    for i, taxon in enumerate(taxa):
        u = lookup_dict[i]
        tree.nodes[u]['label'] = taxon

def all_binary_trees_by_taxa(taxa):
    trees = all_binary_trees_by_num_taxa(len(taxa))
    for tree in trees:
        newtree = tree.copy()
        relabel(newtree, tree.lookup_label_dict, taxa)
        yield newtree
# End of generator usign cached trees

# Sequential generator

def tree_with_single_taxon(taxon):
    t = PhyloTree()
    t.add_node('_1', label=taxon)
    return t

def all_trees(taxa, binary=True):
    """
    Returns a sequential generator for all trees with leaves labeled by `taxa`.

    If `binary` is True, then only binary trees will be generated; otherwise trees will
    have internal nodes with arbitrary out-degree.
    """
    n = len(taxa)
    if n == 1:
        yield tree_with_single_taxon(taxa[0])
        return
    taxon = taxa[-1]
    parent_taxa = taxa[0:-1]
    parent_generator = all_trees(parent_taxa, binary=binary)
    for parent in parent_generator:
        # P&H
        for u in parent.nodes:
            newtree = push_and_hang(parent, u, taxon)
            yield newtree
        if not binary:
            for u in parent.interior_nodes:
                newtree = hold_and_hang(parent, u, taxon)
                yield newtree


@lru_cache(maxsize=None)
def _number_of_trees_bin_global(n):
    """
    Gives the number of binary phylogenetic trees on `n` taxa.
    """
    #return number_of_trees_bin_partial(n, n, 2 * n - 1)
    return product(range(2*n-3,1,-2))


def _random_tree_bin_global(taxa):
    """
    Returns a random binary phylogenetic tree over 'taxa'.
    """

    n = len(taxa)
    if n == 1:
        return tree_with_single_taxon(taxa[0])
    parent = _random_tree_bin_global(taxa[0:-1])
    u = random.choice(list(parent.nodes()))
    newtree = push_and_hang(parent, u, taxa[-1])
    return newtree


# Number of and random trees: not necessarily binary

@lru_cache(maxsize=None)
def _number_of_trees_nobin_partial(n, N):
    """
    Gives the number of multifurcating phylogenetic trees on n taxa and N nodes.
    """
    if (N < n) or (n < 0) or (N >= 2 * n):
        return 0
    if n == 1:
        return 1
    return (N - 2) * _number_of_trees_nobin_partial(n - 1, N - 2) + \
           (N - n) * _number_of_trees_nobin_partial(n - 1, N - 1)


@lru_cache(maxsize=None)
def _number_of_trees_nobin_global(n):
    """
    Gives the number of multifurcating phylogenetic trees on n taxa.
    """
    if n == 1:
        return 1
    return sum([_number_of_trees_nobin_partial(n, N) for N in range(n + 1, 2 * n)])


def _random_tree_nobin_partial(taxa, N):
    """
    Returns a (uniformly) random multifurcating phylogenetic tree over 'taxa' with 'N' nodes.
    """
    n = len(taxa)
    if (N < n) or (n < 0) or (N >= 2 * n):
        return None
    if n == 1:
        return tree_with_single_taxon(taxa[0])
    choices = {
        (push_and_hang, N - 2, 'nodes'): (N - 2) * _number_of_trees_nobin_partial(n - 1, N - 2),
        (hold_and_hang, N - 1, 'interior_nodes'):
            (N - n) * _number_of_trees_nobin_partial(n - 1, N - 1)
    }
    (operation, Np, candidates_method) = random_weighted(choices)
    parent = _random_tree_nobin_partial(taxa[0:-1], Np)
    candidates = list(getattr(parent, candidates_method))
    u = random.choice(candidates)
    newtree = operation(parent, u, taxa[-1])
    return newtree


def _random_tree_nobin_global(taxa):
    """
    Returns a (uniformly) random multifurcating phylogenetic tree over 'taxa'.
    """
    n = len(taxa)
    if n == 1:
        return tree_with_single_taxon(taxa[0])
    choices = {
        N: _number_of_trees_nobin_partial(n, N) for N in range(n + 1, 2 * n)
        }
    N = random_weighted(choices)
    return _random_tree_nobin_partial(taxa, N)


def number_of_trees(n, binary=True):
    """
    Returns the number of phylogenetic trees with `n` labelled leaves.
    If `binary` is True, only  binary trees are counted.
    """
    if binary:
        return _number_of_trees_bin_global(n)
    else:
        return _number_of_trees_nobin_global(n)


def random_uniform_tree_generator(taxa, binary=True):
    """
    Returns a random uniform generator of trees with leaves labeled by `taxa`.

    If `binary` is True, then only binary trees will be generated; otherwise trees will
    have internal nodes with arbitrary out-degree.
    """
    f = _random_tree_bin_global if binary else _random_tree_nobin_global
    while True:
        yield f(taxa)

def random_yule_tree_generator(taxa):
    """
    Returns a random Yule generator of binary trees with leaves labeled by `taxa`.
    """
    ntaxa = len(taxa)
    while True:
        if ntaxa == 0:
            yield None
        if ntaxa == 1:
            yield tree_with_single_taxon(taxa[0])
        tree = PhyloTree()
        u = tree.new_node()
        v = tree.new_node()
        w = tree.new_node()
        tree.add_edges_from([(u,v),(u,w)])
        for label in range(ntaxa - 2):
            u = choice(tree.leaves)
            v = tree.new_node()
            w = tree.new_node()
            tree.add_edges_from([(u, v), (u, w)])
            tree.clear_cache()
        leaves = tree.leaves()
        shuffle(taxa)
        for (leaf, taxon) in zip(leaves, taxa):
            tree.nodes[leaf]['label'] = taxa
        yield tree


if __name__ == "__main__":
    tg = all_trees(['1', '2', '3'])
    while tg:
        print(tg.next())
