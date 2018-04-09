from .classes import PhyloTree
from .operations import push_and_hang, hold_and_hang, push_and_label, hold_and_label
from .utils import random_weighted
from .memoize import memoize_function
from random import randint, shuffle
# def memoize_function(f): return f # use it to document memoized functions (sphinx bug?)


# import numpy
import random


# Sequential generator

def all_trees(taxa, binary=False, nested_taxa=True):
    """
    Returns a generator for trees with given taxa. 
    If nested_taxa is True, then trees with internal nodes labeled 
    will also be produced.
    If binary is True, then only binary trees will be generated; otherwise trees will 
    have internal nodes with arbitrary out-degree. 
    
    EXAMPLE::
    
        >>> generator = all_trees(['A', 'B', 'C'], binary=True)
        >>> generator.next().eNewick()
        ... '((C,B),A);'
        >>> generator.next().eNewick()
        ... '(C,(B,A));'
        >>> tmp = 2
        >>> for tree in generator: tmp += 1
        >>> tmp
        ... 21
    
    EXAMPLE::
    
        >>> for tree in all_trees(['A', 'B', 'C', binary=True, nested_taxa=False):
        >>>     print(tree.eNewick())
        ... ((C,B),A);
        ... (C,(B,A));
        ... ((C,A),B);
        
    EXAMPLE::
    
        >>> for tree in all_trees(['A', 'B', 'C', binary=False, nested_taxa=False):
        >>>     print(tree.eNewick())
        ... ((C,B),A);
        ... (C,(B,A));
        ... ((C,A),B);
        ... (A,B,C);
        
    """
    n = len(taxa)
    if n == 1:
        yield PhyloTree(eNewick=('%s;' % taxa[0]))
        return
    taxon = taxa[-1]
    parent_taxa = taxa[0:-1]
    parent_generator = all_trees(parent_taxa,
                                 binary=binary, nested_taxa=nested_taxa)
    for parent in parent_generator:
        for u in parent.nodes():
            newtree = push_and_hang(parent, u, taxon)
            yield newtree
        #        if not binary:
        for u in parent.nodes():
            if parent.is_leaf(u):
                if nested_taxa:
                    newtree = hold_and_hang(parent, u, taxon)
                    yield newtree
            else:
                if (not binary) or (parent.out_degree(u) == 1):
                    newtree = hold_and_hang(parent, u, taxon)
                    yield newtree
        if nested_taxa:
            for u in parent.nodes():
                newtree = push_and_label(parent, u, taxon)
                yield newtree
            for u in parent.nodes():
                newtree = hold_and_label(parent, u, taxon)
                if newtree:
                    yield newtree


# Number of and random trees: binary without nested taxa

@memoize_function
def number_of_trees_bin_nont_partial(n, l, N):
    """
    Gives the number of phylogenetic trees on n taxa with l leaves and N nodes.
    Assume binary trees without nested taxa.
    """
    if (l != n) or (N != 2 * n - 1) or (n < 0):
        return 0
    if n == 1:
        return 1
    return (N - 2) * number_of_trees_bin_nont_partial(n - 1, l - 1, N - 2)


@memoize_function
def number_of_trees_bin_nont_global(n):
    """
    Gives the number of phylogenetic trees on n taxa.
    Assume binary trees and without nested taxa.
    """
    return number_of_trees_bin_nont_partial(n, n, 2 * n - 1)


def random_tree_bin_nont_global(taxa, id_offset=0):
    """
    Returns a random binary tree without nested taxa.
    
    EXAMPLE::
    
        >>> network = random_tree_bin_nont_global(['A', 'B', 'C', 'D'])
        >>> network.nested_label(network.roots()[0])
        ... '{D,{C,{A,B}}}' # random 
        >>> network = random_tree_bin_nont_global(['A', 'B', 'C', 'D'])
        >>> network.nested_label(network.roots()[0])
        ... '{A,{C,{B,D}}}' # random 
        
    """

    n = len(taxa)
    if n == 1:
        return PhyloTree(eNewick=('%s;' % taxa[0]), id_offset=id_offset)
    parent = random_tree_bin_nont_global(taxa[0:-1], id_offset)
    u = random.choice(parent.nodes())
    newtree = push_and_hang(parent, u, taxa[-1])
    return newtree


# Number of and random trees: not necessarily binary without nested taxa

@memoize_function
def number_of_trees_nobin_nont_partial(n, l, N):
    """
    Gives the number of phylogenetic trees on n taxa with l leaves and N nodes.
    Assume not necessarily binary trees and without nested taxa.
    """
    if (l != n) or (N < n) or (n < 0) or (N >= 2 * n):
        return 0
    if n == 1:
        return 1
    return (N - 2) * number_of_trees_nobin_nont_partial(n - 1, l - 1, N - 2) + \
           (N - n) * number_of_trees_nobin_nont_partial(n - 1, l - 1, N - 1)


@memoize_function
def number_of_trees_nobin_nont_global(n):
    """
    Gives the number of phylogenetic trees on n taxa.
    Assume not necessarily binary trees and without nested taxa.
    """
    if n == 1:
        return 1
    return sum([number_of_trees_nobin_nont_partial(n, n, N) for N in range(n + 1, 2 * n)])


def random_tree_nobin_nont_partial(taxa, l, N, id_offset=0):
    n = len(taxa)
    if (l != n) or (N < n) or (n < 0) or (N >= 2 * n):
        return None
    if n == 1:
        return PhyloTree(eNewick=('%s;' % taxa[0]), id_offset=id_offset)
    choices = {
        (push_and_hang, l - 1, N - 2, 'nodes'): (N - 2) * number_of_trees_nobin_nont_partial(n - 1, l - 1, N - 2),
        (hold_and_hang, l - 1, N - 1, 'interior_nodes'): (N - n) * number_of_trees_nobin_nont_partial(n - 1, l - 1,
                                                                                                      N - 1)
    }
    (operation, lp, Np, candidates_method) = random_weighted(choices)
    parent = random_tree_nobin_nont_partial(taxa[0:-1], lp, Np, id_offset)
    candidates = getattr(parent, candidates_method)()
    u = random.choice(candidates)
    newtree = operation(parent, u, taxa[-1])
    return newtree


def random_tree_nobin_nont_global(taxa, id_offset=0):
    """
    Returns a random tree without nested taxa.
    
    EXAMPLE::
    
        >>> network = random_tree_nobin_nont_global(['A', 'B', 'C', 'D'])
        >>> network.nested_label(network.roots()[0])
        ... '{A,B,{C,D}}' # random, not binary
        >>> network = random_tree_nobin_nont_global(['A', 'B', 'C', 'D'])
        >>> network.nested_label(network.roots()[0])
        ... '{D,{A,{B,C}}}' # random
        
    """
    n = len(taxa)
    if n == 1:
        return PhyloTree(eNewick=('%s;' % taxa[0]), id_offset=id_offset)
    choices = {
        N: number_of_trees_nobin_nont_partial(n, n, N) for N in range(n + 1, 2 * n)
        }
    N = random_weighted(choices)
    l = n
    return random_tree_nobin_nont_partial(taxa, l, N, id_offset)


# Number of and random trees: binary with nested taxa

@memoize_function
def number_of_trees_bin_nt_partial(n, l, N, e):
    """
    Gives the number of phylogenetic trees on n taxa with l leaves, N nodes, e of them being elementary.
    Assume binary trees with nested taxa.
    """
    # print("n=%d, l=%d, N=%d, e=%d" % (n,l,N,e))
    if (l <= 0) or (l > n) or (n < 0) or (N < n) or (l + e > N) or (l + e > n) or (e < 0):
        # h "n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,0)
        return 0
    if n == 1:
        # print("n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,1))
        if (l == 1) and (N == 1) and (e == 0):
            return 1
        else:
            return 0
    count = (
        (N - 2) * number_of_trees_bin_nt_partial(n - 1, l - 1, N - 2, e) +  # P&H
        (N - 1) * number_of_trees_bin_nt_partial(n - 1, l, N - 1, e - 1) +  # P&L
        (N - n + 1) * number_of_trees_bin_nt_partial(n - 1, l, N, e) +  # H&L
        l * number_of_trees_bin_nt_partial(n - 1, l, N - 1, e - 1) +  # H&H leaf
        (e + 1) * number_of_trees_bin_nt_partial(n - 1, l - 1, N - 1, e + 1)  # H&H elem
    )
    # print("n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,count))
    return count


@memoize_function
def number_of_trees_bin_nt_global(n):
    if n == 1:
        return 1
    return sum([number_of_trees_bin_nt_partial(n, l, N, e) for l in range(1, n + 1)
                for N in range(n, 2 * n)
                for e in range(0, n)])


def random_tree_bin_nt_partial(taxa, l, N, e, id_offset):
    n = len(taxa)
    if (l <= 0) or (l > n) or (n < 0) or (N < n) or (l + e > N) or (l + e > n) or (e < 0):
        return None
    if n == 1:
        if (l == 1) and (N == 1) and (e == 0):
            return PhyloTree(eNewick=('%s;' % taxa[0]), id_offset=id_offset)
        else:
            return None
    choices = {
        (push_and_hang, l - 1, N - 2, e, 'nodes'): (N - 2) * number_of_trees_bin_nt_partial(n - 1, l - 1, N - 2, e),
        (push_and_label, l, N - 1, e - 1, 'nodes'): (N - 1) * number_of_trees_bin_nt_partial(n - 1, l, N - 1, e - 1),
        (hold_and_label, l, N, e, 'unlabelled_nodes'): (N - n + 1) * number_of_trees_bin_nt_partial(n - 1, l, N, e),
        (hold_and_hang, l, N - 1, e - 1, 'leaves'): l * number_of_trees_bin_nt_partial(n - 1, l, N - 1, e - 1),
        (hold_and_hang, l - 1, N - 1, e + 1, 'elementary_nodes'): (e + 1) * number_of_trees_bin_nt_partial(n - 1, l - 1,
                                                                                                           N - 1, e + 1)
    }
    (operation, lp, Np, ep, candidates_method) = random_weighted(choices)
    parent = random_tree_bin_nt_partial(taxa[0:-1], lp, Np, ep, id_offset)
    candidates = getattr(parent, candidates_method)()
    u = random.choice(candidates)
    newtree = operation(parent, u, taxa[-1])
    return newtree


def random_tree_bin_nt_global(taxa, id_offset=0):
    """
    Returns a random binary tree with nested taxa.
    
    EXAMPLE::
    
        >>> network = random_tree_bin_nt_global(['A', 'B', 'C', 'D'])
        >>> network.eNewick()
        ... '((D,A),B)C;' # random
        >>> network = random_tree_bin_nt_global(['A', 'B', 'C', 'D'])
        >>> network.eNewick()
        ... '((D,A),(B)C);' # random
        
    """

    n = len(taxa)
    if n == 1:
        return PhyloTree(eNewick=('%s;' % taxa[0]), id_offset=id_offset)
    choices = {
        (l, N, e): number_of_trees_bin_nt_partial(n, l, N, e) for l in range(1, n + 1)
        for N in range(n, 2 * n)
        for e in range(0, n)
        }
    (l, N, e) = random_weighted(choices)
    return random_tree_bin_nt_partial(taxa, l, N, e, id_offset)


# Number of and random trees: not necessarily binary with nested taxa

@memoize_function
def number_of_trees_nobin_nt_partial(n, l, N):
    """
    Gives the number of phylogenetic trees on n taxa with l leaves, N nodes, e of them being elementary.
    Assume binary trees with nested taxa.
    """
    # print("n=%d, l=%d, N=%d, e=%d" % (n,l,N,e))
    if (l <= 0) or (l > n) or (n < 0) or (N < n):
        # print("n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,0))
        return 0
    if n == 1:
        # print("n=%d, l=%d, N=%d, --> %d" % (n,l,N,1))
        if (l == 1) and (N == 1):
            return 1
        else:
            return 0
    count = (
        (N - 2) * number_of_trees_nobin_nt_partial(n - 1, l - 1, N - 2) +  # P&H
        (N - l) * number_of_trees_nobin_nt_partial(n - 1, l - 1, N - 1) +  # H&H interior
        l * number_of_trees_nobin_nt_partial(n - 1, l, N - 1) +  # H&H leaf
        (N - 1) * number_of_trees_nobin_nt_partial(n - 1, l, N - 1) +  # P&L
        (N - n + 1) * number_of_trees_nobin_nt_partial(n - 1, l, N)  # H&L
    )
    # print("n=%d, l=%d, N=%d, --> %d" % (n,l,N,count))
    return count


@memoize_function
def number_of_trees_nobin_nt_global(n):
    if n == 1:
        return 1
    return sum([number_of_trees_nobin_nt_partial(n, l, N) for l in range(1, n + 1)
                for N in range(n, 2 * n)])


def random_tree_nobin_nt_partial(taxa, l, N, id_offset):
    n = len(taxa)
    if (l <= 0) or (l > n) or (n < 0) or (N < n):
        return None
    if n == 1:
        if (l == 1) and (N == 1):
            return PhyloTree(eNewick=('%s;' % taxa[0]), id_offset=id_offset)
        else:
            return None
    choices = {
        (push_and_hang, l - 1, N - 2, 'nodes'): (N - 2) * number_of_trees_nobin_nt_partial(n - 1, l - 1, N - 2),
        (hold_and_hang, l - 1, N - 1, 'interior_nodes'): (N - l) * number_of_trees_nobin_nt_partial(n - 1, l - 1,
                                                                                                    N - 1),
        (hold_and_hang, l, N - 1, 'leaves'): l * number_of_trees_nobin_nt_partial(n - 1, l, N - 1),
        (push_and_label, l, N - 1, 'nodes'): (N - 1) * number_of_trees_nobin_nt_partial(n - 1, l, N - 1),
        (hold_and_label, l, N, 'unlabelled_nodes'): (N - n + 1) * number_of_trees_nobin_nt_partial(n - 1, l, N)
    }
    # print(choices)
    (operation, lp, Np, candidates_method) = random_weighted(choices)
    parent = random_tree_nobin_nt_partial(taxa[0:-1], lp, Np, id_offset)
    # print(parent)
    candidates = getattr(parent, candidates_method)()
    u = random.choice(candidates)
    newtree = operation(parent, u, taxa[-1])
    return newtree


def random_tree_nobin_nt_global(taxa, id_offset):
    """
    Returns a random tree with nested taxa.
    
    EXAMPLE::
    
        >>> network = random_tree_nobin_nt_global(['A', 'B', 'C', 'D'])
        >>> network.eNewick()
        ... '((D,A),C,B);' # random, non-binary
        >>> network = random_tree_nobin_nt_global(['A', 'B', 'C', 'D'])
        >>> network.eNewick()
        ... '(((D,C),B),A);' # random
        
    """

    n = len(taxa)
    if n == 1:
        return PhyloTree(eNewick=('%s;' % taxa[0]), id_offset=id_offset)
    choices = {
        (l, N): number_of_trees_nobin_nt_partial(n, l, N) for l in range(1, n + 1)
        for N in range(n, 2 * n)
        }
    # print(choices)
    (l, N) = random_weighted(choices)
    # print(l,N)
    return random_tree_nobin_nt_partial(taxa, l, N, id_offset)


def random_yule_tree(taxa, id_offset=0):
    ntaxa = len(taxa)
    if ntaxa == 0:
        return None
    if ntaxa == 1:
        return PhyloTree(eNewick=taxa[0] + ";")
    network = PhyloTree()
    network.add_node('_' + str(id_offset))
    network.add_node('_' + str(id_offset + 1))
    network.add_node('_' + str(id_offset + 2))
    network.add_edge('_' + str(id_offset), '_' + str(id_offset + 1))
    network.add_edge('_' + str(id_offset), '_' + str(id_offset + 2))
    id_offset += 3
    for label in range(ntaxa - 2):
        nn = len(network.leaves())
        parent = network.leaves()[randint(0, nn - 1)]

        newnode = '_' + str(id_offset)
        network.add_node(newnode)
        network.add_edge(parent, newnode)
        newnode = '_' + str(id_offset + 1)
        network.add_node(newnode)
        network.add_edge(parent, newnode)

        id_offset += 2
        network.cache = {}
    leaves = network.leaves()
    shuffle(taxa)
    for i in range(ntaxa):
        network.set_label(leaves[i], taxa[i])
    return network


def random_tree_generator(taxa, binary=False, nested_taxa=True, yule=False, id_offset=0):
    """
    Returns generator of random trees.
    If nested_taxa = True, then trees with internal labels will also be produced.
    If binary = True, then only binary trees will be produced.
    
    If yule = True, instead of generating trees with uniform probability we will 
    use the yule model. Remember that yule model is only valid for non-nested 
    binary trees.
    
    EXAMPLE::
    
        >>> generator = random_tree_generator(['a', 'b', 'c', 'd', 'e', 'f'], binary=False)
        >>> for i in range(3):
        >>>     print(generator.next().eNewick())
        ... ((a)d,(e,b,(c,f)));
        ... (((((d)c)e)f,a))b;
        ... ((e,d,(c)a),b,f);
        
    EXAMPLE::
    
        >>> generator = random_tree_generator(['a', 'b', 'c', 'd', 'e', 'f'], binary=True)
        >>> for i in range(3):
        >>>     print(generator.next().eNewick())
        ... '((d,b),(c,((e)f)a));'
        ... '(c,(((a,e),d),f)b);'
        ... '(((c,f),((b)e)d))a;'
        
    EXAMPLE::
    
        >>> generator = random_tree_generator(['a', 'b', 'c', 'd', 'e', 'f'], yule=True)
        >>> for i in range(3):
        >>>     print(generator.next().eNewick())
        ... ((b,(d,f)),((e,c),a));
        ... (((c,e),b),((a,f),d));
        ... ((e,(f,(a,b))),(d,c));
        
    """

    if binary and not yule:
        if nested_taxa:
            f = random_tree_bin_nt_global
        else:
            f = random_tree_bin_nont_global
    elif yule:
        f = random_yule_tree
    else:
        if nested_taxa:
            f = random_tree_nobin_nt_global
        else:
            f = random_tree_nobin_nont_global
    while True:
        yield f(taxa, id_offset)


if __name__ == "__main__":
    tg = all_trees(['1', '2', '3'])
    while tg:
        print(tg.next())
