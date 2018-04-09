"""
The phylonetwork is the main module of the project PhyloNetworks.

It includes classes to represent phylogenetic networks and trees with or
without nested taxa, as well as functions for both the successive and
random generation of all phylogenetic trees over a set of taxa.
"""

from .classes import *
from .distances import *
from .memoize import *
from .permutations import *
from .eNewick import *
from .generators import *
from .exceptions import *

def eNewick_reader(filename,ignore_prefix=None):
    f = open(filename,'r')
    while True:
        line=f.readline().strip()
        if not line: return
        net=PhyloNetwork(eNewick=line,ignore_prefix=ignore_prefix)
        yield net


    
if __name__ == '__main__':
    net1=PhyloNetwork(eNewick='(1,2,3);')
    net2=PhyloNetwork(eNewick='(1,2)3;')
    print(mu_distance(net1,net2))
    
