from _thread import RLock
from random import randint


def random_weighted(d):
    """Given ``d`` a dictionary, returns a key in it with probability equals to its value 
    (normalized over the sum of all values)
    """
    sum_of_values = sum(d.values())
    r = randint(0, sum_of_values-1)
    for item in d:
        if r < d[item]:
            return item
        r = r-d[item]


def product(iterable):
    result = 1
    for item in iterable:
        result *= item
    return result


_NOT_FOUND = object()


class clearable_cached_property:
    """
    Adapted from standard library in order to keep track of cached attributes and allow for its clearance.
    """
    def __init__(self, func):
        self.func = func
        self.attrname = None
        self.__doc__ = func.__doc__
        self.lock = RLock()

    def __set_name__(self, owner, name):
        if self.attrname is None:
            self.attrname = name
        elif name != self.attrname:
            raise TypeError(
                "Cannot assign the same cached_property to two different names "
                f"({self.attrname!r} and {name!r})."
            )

    def __get__(self, instance, owner=None):
        if instance is None:
            return self
        if self.attrname is None:
            raise TypeError(
                "Cannot use cached_property instance without calling __set_name__ on it.")
        try:
            cache = instance.__dict__
        except AttributeError:  # not all objects have __dict__ (e.g. class defines slots)
            msg = (
                f"No '__dict__' attribute on {type(instance).__name__!r} "
                f"instance to cache {self.attrname!r} property."
            )
            raise TypeError(msg) from None
        val = cache.get(self.attrname, _NOT_FOUND)
        if val is _NOT_FOUND:
            with self.lock:
                # check if another thread filled cache while we awaited lock
                val = cache.get(self.attrname, _NOT_FOUND)
                if val is _NOT_FOUND:
                    val = self.func(instance)
                    try:
                        cache[self.attrname] = val
                        try:
                            instance.__cached_attributes__.append(self.attrname)
                        except AttributeError:
                            instance.__cached_attributes__ = [self.attrname]
                    except TypeError:
                        msg = (
                            f"The '__dict__' attribute on {type(instance).__name__!r} instance "
                            f"does not support item assignment for caching {self.attrname!r} property."
                        )
                        raise TypeError(msg) from None
        return val


def clear_cache(obj):
    try:
        cached_attributes = obj.__cached_attributes__
    except AttributeError:
        return
    for attr in cached_attributes:
        delattr(obj, attr)
    obj.__cached_attributes__ = []