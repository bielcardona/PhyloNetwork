import copy


class Permutation():
    """
	Class for permutations
	"""

    def __init__(self, mapping=None, ordering=None):
        self.mapping = {}

        if not ordering == None:
            keys = ordering[:]
            keys.sort()
            for i in range(len(keys)):
                self.mapping[keys[i]] = ordering[i]

        if not mapping == None:
            self.mapping = copy.copy(mapping)

    def __str__(self):
        return str(self.mapping)

    def __getitem__(self, key):
        try:
            return self.mapping[key]
        except:
            return key

    def __setitem__(self, key, value):
        self.mapping[key] = value

    def __mul__(self, perm2):
        result = Permutation()
        keys = set(perm2.mapping.keys())
        keys.update(list(self.mapping.keys()))
        for key in keys:
            try:
                result[key] = self[perm2[key]]
            except:
                result[key] = perm2[key]
        return result

    def inverse(self):
        result = Permutation()
        for key in self.mapping:
            result[self[key]] = key
        return result

    def __pow__(self, n):
        if n < 0:
            inv = self.inverse()
            result = copy.copy(inv)
            for i in range(-n - 1):
                result = result * inv
            return result
        elif n > 0:
            result = copy.copy(self)
            for i in range(n - 1):
                result = result * self
            return result
        else:
            return Permutation()

    def identity(self, keys):
        for key in keys:
            self[key] = key

    def cycles(self):
        keys = set(self.mapping.keys())
        cyclesfound = []
        while keys:
            key = keys.pop()
            cycle = [key]
            actual = self[key]
            while not actual == key:
                cycle.append(actual)
                keys.remove(actual)
                actual = self[actual]
            cyclesfound.append(cycle)
        return cyclesfound
