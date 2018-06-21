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

def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K
