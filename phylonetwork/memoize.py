import functools


class memoize_method(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated. The cached values are stored in the cache property of the
    instance whose method is being decorated.
    """

    def __init__(self, func):
        # print("Init")
        self.func = func

    def __call__(self, *args):
        # print("Call")
        if not self.func in self.cache:
            self.cache[self.func] = {}
        try:
            return self.cache[self.func][args]
        except KeyError:
            value = self.func(*args)
            self.cache[self.func][args] = value
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """Support instance methods."""
        # print("Get", obj, objtype)
        fn = functools.partial(self.__call__, obj)
        fn.__doc__ = self.func.__doc__
        fn.__name__ = self.func.__name__
        try:
            self.cache = obj.cache
        except:
            obj.cache = {}
            self.cache = obj.cache
        # print(self.cache)
        return fn


class memoize_function(object):
    """Decorator that caches a function's return value each time it is called.
    If called later with the same arguments, the cached value is returned, and
    not re-evaluated.
    """

    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        try:
            return self.cache[args]
        except KeyError:
            value = self.func(*args)
            self.cache[args] = value
            return value
        except TypeError:
            # uncachable -- for instance, passing a list as an argument.
            # Better to not cache than to blow up entirely.
            return self.func(*args)

    def __repr__(self):
        """Return the function's docstring."""
        return self.func.__doc__

    def __get__(self, obj, objtype):
        """Support instance methods."""
        return functools.partial(self.__call__, obj)


if __name__ == "__main__":
    class MyClass(object):
        def __init__(self, data):
            self.data = data

        def update(self, data):
            self.data = data
            self.cache = {}

        @memoize_method
        def func1(self, x):
            print("Computing func1")
            return "I am func1 of %s. Data is %s. x is %s\n" % (self, self.data, x)

        @memoize_method
        def func2(self, x):
            print( "Computing func2")
            return "I am func2 of %s. Data is %s. x is %s\n" % (self, self.data, x)

        def func3(self, x):
            print( "Computing func3")
            return "I am func3 of %s. Data is %s. x is %s\n" % (self, self.data, x)


    mc1 = MyClass("data1")
    mc2 = MyClass("data2")
    mc3 = MyClass("data3")

    print( mc1.func1(1))
    print( mc1.func1(1))
    print( mc1.func2(1))
    print( mc1.func2(1))
    print( mc1.func3(1))
    print( mc1.func3(1))

    print( mc2.func1(1))
    print( mc2.func1(1))
    print( mc2.func2(1))
    print( mc2.func2(1))
    print( mc2.func3(1))
    print( mc2.func3(1))

    print( "Update mc1\n")
    mc1.update("data1new")

    print( mc1.func1(1))
    print( mc1.func2(1))
    print( mc1.func3(1))
    print( mc2.func1(1))
    print( mc2.func2(1))
    print( mc2.func3(1))
