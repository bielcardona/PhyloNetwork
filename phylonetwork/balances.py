"""
Implementation of different balance indices for phylogenetic trees or networks.
"""
import phylonetwork


def sackin(t):
    """Computes the Sackin (un)balance index of `t`"""
    return sum([t.depths[u] for u in t.leaves])

def extended_sackin(net:phylonetwork.PhylogeneticNetwork):
    """Computes the extended Sackin (un)balance index of `net`"""
    factor = sum(net.mu(net.root))
    weights = [sum(net.mu(v)) for v in net.tree_nodes]
    return sum(weights)/factor


def entropy(t):
    """Computes the entropy balance index of `t`"""
    return sum([t.depths[u] * 2**(-t.depths[u]) for u in t.leaves])
