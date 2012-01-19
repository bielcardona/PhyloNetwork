'''
Created on Jan 16, 2012

@author: cardona
'''

from .classes import PhyloNetwork
from .distances import mu_distance
from .operations import push_and_hang,hold_and_hang,push_and_label,hold_and_label

#import numpy
import random

def Tree_generator(taxa,binary=False,nested_taxa=True):
    """Returns a generator for trees with given taxa. 
    If nested_taxa is True, then trees with internal nodes labeled 
    will also be produced.
    If binary is True, then only binary trees will be generated; otherwise trees will 
    have internal nodes with arbitrary out-degree. 
    """
    n=len(taxa)
    if n==1:
        yield PhyloNetwork(eNewick=('%s;' % taxa[0]))
        return
    taxon=taxa[-1]
    parent_taxa=taxa[0:-1]
    parent_generator=Tree_generator(parent_taxa,
                                    binary=binary,nested_taxa=nested_taxa)
    for parent in parent_generator:
        for u in parent.nodes():
            newtree=push_and_hang(parent,u,taxon)
            yield newtree
#        if not binary:
        for u in parent.nodes():
            if parent.is_leaf(u):
                if nested_taxa:
                    newtree=hold_and_hang(parent,u,taxon)
                    yield newtree
            else:
                if (not binary) or (parent.out_degree(u)==1):
                    newtree=hold_and_hang(parent,u,taxon)
                    yield newtree
        if nested_taxa:
            for u in parent.nodes():
                newtree=push_and_label(parent,u,taxon)
                yield newtree
            for u in parent.nodes():
                newtree=hold_and_label(parent,u,taxon)
                if newtree:
                    yield newtree

def number_of_trees_bin_nont_partial(n,l,N):
    """Gives the number of phylogenetic trees on n taxa with l leaves and N nodes.
    Assume binary trees without nested taxa.
    """
    if (l != n) or (N != 2*n-1) or (n < 0):
        return 0
    if n == 1:
        return 1
    return (N-2)*number_of_trees_bin_nont_partial(n-1,l-1,N-2)

def number_of_trees_bin_nont_global(n):
    """Gives the number of phylogenetic trees on n taxa.
    Assume binary trees and without nested taxa.
    """
    return number_of_trees_bin_nont_partial(n,n,2*n-1)

def number_of_trees_nobin_nont_partial(n,l,N):
    """Gives the number of phylogenetic trees on n taxa with l leaves and N nodes.
    Assume not necessarily binary trees and without nested taxa.
    """
    if (l != n) or (N < n) or (n < 0) or (N >= 2*n):
        return 0
    if n==1:
        return 1
    return (N-2)*number_of_trees_nobin_nont_partial(n-1,l-1,N-2) + \
           (N-n)*number_of_trees_nobin_nont_partial(n-1,l-1,N-1)

def number_of_trees_nobin_nont_global(n):
    """Gives the number of phylogenetic trees on n taxa.
    Assume not necessarily binary trees and without nested taxa.
    """
    if n == 1:
        return 1
    return sum([number_of_trees_nobin_nont_partial(n,n,N) for N in range(n+1,2*n)])

def number_of_trees_bin_nt_partial(n,l,N,e,log=False):
    """Gives the number of phylogenetic trees on n taxa with l leaves, N nodes, e of them being elementary.
    Assume binary trees with nested taxa.
    """
    #print "n=%d, l=%d, N=%d, e=%d" % (n,l,N,e)
    if (l <= 0) or (l > n) or (n < 0) or (N < n) or (l+e > N) or (l+e > n) or (e < 0):
        #print "n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,0)
        return 0
    if n==1:
        #print "n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,1)
        if (l == 1) and (N == 1) and (e == 0):
            return 1
        else:
            return 0
    count = (
        (N-2)*number_of_trees_bin_nt_partial(n-1,l-1,N-2,e) + # P&H
        (N-1)*number_of_trees_bin_nt_partial(n-1,l,N-1,e-1) + # P&L
        (N-n+1)*number_of_trees_bin_nt_partial(n-1,l,N,e) +   # H&L
        (l)*number_of_trees_bin_nt_partial(n-1,l,N-1,e-1) +   # H&H leaf
        (e+1)*number_of_trees_bin_nt_partial(n-1,l-1,N-1,e+1) # H&H elem
    )
    #print "n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,count)
    return count

def number_of_trees_bin_nt_global(n):
    if n==1:
        return 1
    return sum([number_of_trees_bin_nt_partial(n,l,N,e) for l in range(1,n+1) \
                                                        for N in range(n,2*n)
                                                        for e in range(0,n)])

def number_of_trees_nobin_nt_partial(n,l,N,log=False):
    """Gives the number of phylogenetic trees on n taxa with l leaves, N nodes, e of them being elementary.
    Assume binary trees with nested taxa.
    """
    #print "n=%d, l=%d, N=%d, e=%d" % (n,l,N,e)
    if (l <= 0) or (l > n) or (n < 0) or (N < n):
        #print "n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,0)
        return 0
    if n==1:
        #print "n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,1)
        if (l == 1) and (N == 1): 
            return 1
        else:
            return 0
    count = (
        (N-2)*number_of_trees_nobin_nt_partial(n-1,l-1,N-2) + # P&H
        (N-1)*number_of_trees_nobin_nt_partial(n-1,l-1,N-1) + # H&H
        (N-1)*number_of_trees_nobin_nt_partial(n-1,l,N-1) +   # P&L
        (N-n+1)*number_of_trees_nobin_nt_partial(n-1,l,N)   # H&L 
    )
    #print "n=%d, l=%d, N=%d, e=%d, --> %d" % (n,l,N,e,count)
    return count

def number_of_trees_nobin_nt_global(n):
    if n==1:
        return 1
    return sum([number_of_trees_nobin_nt_partial(n,l,N) for l in range(1,n+1) \
                                                        for N in range(n,2*n)])

_not_byint={}
_not_byint[(1,0)]=1
_not_global={}

def number_of_trees_numint(n,k):
    """Returns the number of phylogenetic trees with n taxa and with k internal nodes
    """
    try:
        return _not_byint[(n,k)]
    except:
        if k<=0 or k>=n:
            _not_byint[(n,k)]=0
        elif n==1 or k==1:
            _not_byint[(n,k)]=1
        else:
            _not_byint[(n,k)]=k * number_of_trees_numint(n-1,k) + \
                (n+k-2)* number_of_trees_numint(n-1,k-1)
    return _not_byint[(n,k)]

def number_of_trees_global(n):
    try:
        return _not_global[n]
    except:
        _not_global[n]=sum([number_of_trees_numint(n,k) for k in range(n)])
    return _not_global[n]

def rtg_internal(taxa,k):
    n=len(taxa)
    if n==1 and k==0:
        return PhyloNetwork(eNewick=('%s;' % taxa[0]))
    nombre=number_of_trees_numint(n,k)
    aleatori=random.randint(1,nombre)
    if aleatori <= (n+k-2)*number_of_trees_numint(n-1,k-1):
        # comes from push_and_hang
        parent=rtg_internal(taxa[0:-1],k-1)
        node=random.choice(parent.nodes())
        newtree=push_and_hang(parent,node,taxa[-1])
    else:
        # comes from hold_and_hang
        parent=rtg_internal(taxa[0:-1],k)
        nodes=[u for u in parent.nodes() if not parent.is_leaf(u)]
        node=random.choice(nodes)
        newtree=hold_and_hang(parent,node,taxa[-1])
    return newtree
                    
def random_tree_generator(taxa,binary=False): # no nested taxa implemented
    n=len(taxa)
    numbers=[number_of_trees_numint(n,k) for k in range(n)]
    total=number_of_trees_global(n)
    while True:
        rnd=random.randint(1,total)
        k=1
        while rnd>numbers[k]:
            rnd -= numbers[k]
            k += 1
        yield rtg_internal(taxa,k)
        
def test_uniform(taxa,num):
    alltrees=list(Tree_generator(taxa,nested_taxa=False))
    numtrees=len(alltrees)
    stat={}
    randomfac=random_tree_generator(taxa)
    for i in range(num):
        tree=randomfac.next()
        trobat=False
        for j in range(numtrees):
            if mu_distance(alltrees[j],tree) == 0:
                trobat=True
                try:
                    stat[j] += 1
                except:
                    stat[j] = 1
                break
        if not trobat:
            print 'error'
            print tree.eNewick()
    return stat

if __name__ == "__main__":
    tg = Tree_generator(['1','2','3'])
    while tg:
        print tg.next()
        
