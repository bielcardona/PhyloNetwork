'''
Created on Dec 24, 2011

@author: cardona
'''

def total_cmp(x,y):
    """
    Return -1 if x<y, 0 if x=y, 1 if x>y, with respect to the product partial
    order.
    """
    n=x.size
    for i in range(n):
        if x[i]<y[i]: return -1
        if x[i]>y[i]: return 1
    return 0
