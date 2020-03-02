from random import randint


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


def product(iterable):
    result = 1
    for item in iterable:
        result *= item
    return result

