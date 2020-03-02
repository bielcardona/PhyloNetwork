from phylonetwork import NetworkShape, PhylogeneticNetwork
import random
from itertools import combinations


def speciate(net: NetworkShape, u):
    result = net.copy()
    v = next(result.predecessors(u))
    w = result.make_elementary_node(v, u)
    up = result.new_node()
    result.add_edge(w, up)
    result.clear_cache()
    return result


def hybridize(net: NetworkShape, u1, u2):
    result = net.copy()
    v1 = next(result.predecessors(u1))
    v2 = next(result.predecessors(u2))
    w1 = result.make_elementary_node(v1, u1)
    w2 = result.make_elementary_node(v2, u2)
    h = result.new_reticulation_node()
    l = result.new_node()
    result.add_edges_from([(w1, h), (w2, h), (h, l)])
    result.clear_cache()
    return result


def random_combination(base, r):
    result = set()
    while len(result) != r:
        result = set()
        for _ in range(r):
            result.add(random.choice(base))
    return tuple(result)


def random_yule_btctc_generator(taxa, p):
    ntaxa = len(taxa)
    while True:
        if ntaxa == 0:
            yield None
        if ntaxa == 1:
            yield PhylogeneticNetwork(eNewick=f"{taxa[0]};")
        if ntaxa == 2:
            yield PhylogeneticNetwork(eNewick=f"({taxa[0]},{taxa[1]});")
        net = PhylogeneticNetwork()
        r = net.new_node()
        u = net.new_node()
        v = net.new_node()
        net.add_edges_from([(r, u), (r, v)])
        for _ in range(ntaxa - 2):
            if random.random() < p:
                # make speciate
                u = random.choice(list(net.leaves))
                net = speciate(net, u)
            else:
                # make hybridize
                u, v = random_combination(list(net.leaves), 2)
                net = hybridize(net, u, v)
        leaves = list(net.leaves)
        random.shuffle(taxa)
        for (leaf, taxon) in zip(leaves, taxa):
            net.nodes[leaf]['label'] = taxon
        yield net




