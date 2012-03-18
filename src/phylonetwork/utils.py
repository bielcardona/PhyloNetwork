from random import randint

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

def random_weighted(d):
    """Given ``d`` a dictionary, returns a key in it with probability equals to its value 
    (normalized over the sum of all values)
    """
    sum_of_values = sum(d.values())
    r = randint(0,sum_of_values-1)
    for item in d:
        if r < d[item]:
            return item
        r = r-d[item]
        