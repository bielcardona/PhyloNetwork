"""
Class for phylogenetic networks... with nested taxa
"""

from classes import *

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
    print net1.mu_distance(net2)
    
