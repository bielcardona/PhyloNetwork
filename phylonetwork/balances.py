def sackin(t):
    return sum([t.depths[u] for u in t.leaves])


def entropy(t):
    return sum([t.depths[u] * 2**(-t.depths[u]) for u in t.leaves])
