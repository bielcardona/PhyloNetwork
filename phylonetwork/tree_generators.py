from . import PhylogeneticTree as PhyloTree
from .tree_operations import push_and_hang, hold_and_hang, push_and_label, hold_and_label
from .utils import random_weighted, product
from random import randint, shuffle, choice
from functools import lru_cache
import random


# Sequential generator

def all_trees(taxa, binary=True):
    """
    Returns a generator for trees with given taxa. 

    If binary is True, then only binary trees will be generated; otherwise trees will
    have internal nodes with arbitrary out-degree.
    """
    n = len(taxa)
    if n == 1:
        yield PhyloTree(Newick=('%s;' % taxa[0]))
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
def number_of_trees_bin_global(n):
    """
    Gives the number of binary phylogenetic trees on 'n' taxa.
    """
    #return number_of_trees_bin_partial(n, n, 2 * n - 1)
    return product(range(2*n-3,1,-2))


def random_tree_bin_global(taxa):
    """
    Returns a random binary phylogenetic tree over 'taxa'.
    """

    n = len(taxa)
    if n == 1:
        return PhyloTree(Newick=('%s;' % taxa[0]))
    parent = random_tree_bin_global(taxa[0:-1])
    u = random.choice(list(parent.nodes()))
    newtree = push_and_hang(parent, u, taxa[-1])
    return newtree


# Number of and random trees: not necessarily binary

@lru_cache(maxsize=None)
def number_of_trees_nobin_partial(n, N):
    """
    Gives the number of multifurcating phylogenetic trees on n taxa and N nodes.
    """
    if (N < n) or (n < 0) or (N >= 2 * n):
        return 0
    if n == 1:
        return 1
    return (N - 2) * number_of_trees_nobin_partial(n - 1, N - 2) + \
           (N - n) * number_of_trees_nobin_partial(n - 1, N - 1)


@lru_cache(maxsize=None)
def number_of_trees_nobin_global(n):
    """
    Gives the number of multifurcating phylogenetic trees on n taxa.
    """
    if n == 1:
        return 1
    return sum([number_of_trees_nobin_partial(n, N) for N in range(n + 1, 2 * n)])


def random_tree_nobin_partial(taxa, N):
    """
    Returns a (uniformly) random multifurcating phylogenetic tree over 'taxa' with 'N' nodes.
    """
    n = len(taxa)
    if (N < n) or (n < 0) or (N >= 2 * n):
        return None
    if n == 1:
        return PhyloTree(Newick=('%s;' % taxa[0]))
    choices = {
        (push_and_hang, N - 2, 'nodes'): (N - 2) * number_of_trees_nobin_partial(n - 1, N - 2),
        (hold_and_hang, N - 1, 'interior_nodes'):
            (N - n) * number_of_trees_nobin_partial(n - 1, N - 1)
    }
    (operation, Np, candidates_method) = random_weighted(choices)
    parent = random_tree_nobin_partial(taxa[0:-1], Np)
    candidates = list(getattr(parent, candidates_method))
    u = random.choice(candidates)
    newtree = operation(parent, u, taxa[-1])
    return newtree


def random_tree_nobin_global(taxa):
    """
    Returns a (uniformly) random multifurcating phylogenetic tree over 'taxa'.
    """
    n = len(taxa)
    if n == 1:
        return PhyloTree(Newick=('%s;' % taxa[0]))
    choices = {
        N: number_of_trees_nobin_partial(n, N) for N in range(n + 1, 2 * n)
        }
    N = random_weighted(choices)
    return random_tree_nobin_partial(taxa, N)


def random_uniform_tree_generator(taxa, binary=True):
    f = random_tree_bin_global if binary else random_tree_nobin_global
    while True:
        yield f(taxa)

def random_yule_tree_generator(taxa):
    ntaxa = len(taxa)
    while True:
        if ntaxa == 0:
            yield None
        if ntaxa == 1:
            yield PhyloTree(Newick=taxa[0] + ";")
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
