"""
Implementation of different balance indices for phylogenetic trees or networks.
"""

def sackin(t):
    """Computes the Sackin (un)balance index of `t`"""
    return sum([t.depths[u] for u in t.leaves])


def entropy(t):
    """Computes the entropy balance index of `t`"""
    return sum([t.depths[u] * 2**(-t.depths[u]) for u in t.leaves])
