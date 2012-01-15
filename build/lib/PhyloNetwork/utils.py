'''
Created on Dec 24, 2011

@author: cardona
'''

def total_cmp(x,y):
    """
    Return -1 if x<y, 0 if x=y, 1 if x>y, with respect to the product partial
    order.
    """
    nx=x.size
    ny=y.size
    if nx < ny: return -1
    if nx > ny: return 1
    for i in range(nx):
        if x[i]<y[i]: return -1
        if x[i]>y[i]: return 1
    return 0
