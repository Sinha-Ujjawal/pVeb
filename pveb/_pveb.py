from math import sqrt, ceil, log2
from bisect import bisect_left, bisect_right

from pyrsistent import pmap
from pyrsistent._compat import Hashable, Mapping


class PVeb:
    class _Node:
        def __init__(self, x, u, threshold, data):
            self.min = self.max = x
            self.min_data = self.max_data = data

            self.u = u
            self.threshold = threshold

            if u > threshold:
                self.cluster = pmap()
                self.summary = None
            else:
                self.map = pmap({x: data})

        def copy(self):
            new = PVeb._Node(self.min, self.u, self.threshold, self.min_data)

            new.min = self.min
            new.min_data = self.min_data

            new.max = self.max
            new.max_data = self.max_data

            if new.u > new.threshold:
                new.cluster = self.cluster
                new.summary = self.summary
            else:
                new.map = self.map

            return new

    def __init__(self, lb=0, ub=(1 << 32) - 1, c=100):
        """Constructor for persistent van-emde-boas
        
        Keyword Arguments:
            lb {int} -- [lb of the van-emde-boas tree] (default: {0})
            ub {int} -- [ub of the van-emde-boas tree] (default: {(1 << 32)-1})
            c {int} -- [
                parameter to define van-emde-boas height threshold, higher the c,
                lower the height, hence slower can be the search. Choose wisely
            ] (default: {100})
        """
        assert (
            isinstance(lb, int)
            and isinstance(ub, int)
            and lb <= ub
            and isinstance(c, int)
            and c > 0
        )
        self._lb = lb
        self._ub = ub
        self._c = c

        self._u = ub - lb + 1
        self._threshold = max(2, ceil(c * log2(log2(self._u))))

        self._root = None
        self._len = 0

    @property
    def lb(self):
        return self._lb

    @property
    def ub(self):
        return self._ub

    @property
    def c(self):
        return self._c

    @property
    def u(self):
        return self._u

    @property
    def threshold(self):
        return self._threshold

    def __len__(self):
        return self._len

    def __iter__(self):
        return self.iter_items()

    def __repr__(self):
        return "pvEB({})".format({k: v for k, v in self})

    def __contains__(self, x):
        return self.contains(x)

    def __getitem__(self, x):
        return self.get(x)

    def _check_x(self, x):
        assert isinstance(x, int) and self._lb <= x <= self._ub

    @staticmethod
    def _hierarchical_coords(x, u):
        return x // u, x % u

    @staticmethod
    def _index_of(i, j, u):
        return i * u + j

    @staticmethod
    def _set(V, x, u, threshold, data):
        if V == None:
            return PVeb._Node(x, u, threshold, data)
        if x == V.min or x == V.max:
            return V
        V = V.copy()
        if u > threshold:
            if x < V.min:
                x, V.min = V.min, x
                data, V.min_data = V.min_data, data
            elif x > V.max:
                V.max, V.max_data = x, data

            new_u = ceil(sqrt(u))
            c, i = PVeb._hierarchical_coords(x, new_u)

            if c not in V.cluster:
                V.summary = PVeb._set(V.summary, c, new_u, threshold, data)
            V.cluster = V.cluster.set(
                c, PVeb._set(V.cluster.get(c, None), i, new_u, threshold, data)
            )
        else:
            if x < V.min:
                V.min, V.min_data = x, data
            elif x > V.max:
                V.max, V.max_data = x, data
            V.map = V.map.set(x, data)
        return V

    @staticmethod
    def _pop(V, x, u, threshold):
        if V:
            V = V.copy()
            if u > threshold:
                new_u = ceil(sqrt(u))
                if x == V.min:
                    if V.summary != None:
                        c = V.summary.min
                        x = V.min = PVeb._index_of(c, V.cluster[c].min, new_u)
                        V.min_data = V.cluster[c].min_data
                    else:
                        return None
                c, i = PVeb._hierarchical_coords(x, new_u)
                t = PVeb._pop(V.cluster[c], i, new_u, threshold)
                if t == None:
                    V.cluster = V.cluster.pop(c)
                    V.summary = PVeb._pop(V.summary, c, new_u, threshold)
                if V.summary == None:
                    V.max = V.min
                    V.max_data = V.min_data
                else:
                    c = V.summary.max
                    V.max = PVeb._index_of(c, V.cluster[c].max, new_u)
                    V.max_data = V.cluster[c].max_data
            else:
                V.map = V.map.pop(x)

        return V

    @staticmethod
    def _get(V, x, u, threshold):
        if V:
            if u > threshold:
                if x == V.min:
                    return True, V.min_data
                elif x == V.max:
                    return True, V.max_data
                else:
                    new_u = ceil(sqrt(u))
                    c, i = PVeb._hierarchical_coords(x, new_u)
                    return PVeb._get(V.cluster.get(c, None), i, new_u, threshold)
            else:
                if x in V.map:
                    return True, V.map[x]
        return False, None

    @staticmethod
    def _successor(V, x, u, threshold):
        if V == None:
            return None
        if u > threshold:
            if x < V.min:
                return V.min, V.min_data
            new_u = ceil(sqrt(u))
            c, i = PVeb._hierarchical_coords(x, new_u)
            if c in V.cluster and i < V.cluster[c].max:
                i, data = PVeb._successor(V.cluster[c], i, new_u, threshold)
                return PVeb._index_of(c, i, new_u), data
            else:
                c = PVeb._successor(V.summary, c, new_u, threshold)
                if c != None:
                    c, data = c
                    return PVeb._index_of(c, V.cluster[c].min, new_u), data
        else:
            if len(V.map):
                keys = sorted(V.map.keys())
                index = bisect_right(keys, x)
                if index < len(V.map):
                    return keys[index], V.map[keys[index]]

    @staticmethod
    def _predecessor(V, x, u, threshold):
        if V == None:
            return None
        if u > threshold:
            if x > V.max:
                return V.max, V.max_data
            new_u = ceil(sqrt(u))
            c, i = PVeb._hierarchical_coords(x, new_u)
            if c in V.cluster and i > V.cluster[c].min:
                i, data = PVeb._predecessor(V.cluster[c], i, new_u, threshold)
                return PVeb._index_of(c, i, new_u), data
            else:
                c = PVeb._predecessor(V.summary, c, new_u, threshold)
                if c != None:
                    c, data = c
                    return PVeb._index_of(c, V.cluster[c].max, new_u), data
        else:
            if len(V.map):
                keys = sorted(V.map.keys())
                index = bisect_left(keys, x)
                if index:
                    index -= 1
                    return keys[index], V.map[keys[index]]

    def set(self, x, data):
        self._check_x(x)
        if x not in self:
            root = PVeb._set(self._root, x - self._lb, self._u, self._threshold, data)
            new = PVeb(lb=self._lb, ub=self._ub, c=self._c)
            new._root = root
            new._len = self._len + 1
            return new
        else:
            return self

    def pop(self, x):
        self._check_x(x)
        if x in self:
            root = PVeb._pop(self._root, x - self._lb, self._u, self._threshold)
            new = PVeb(lb=self._lb, ub=self._ub, c=self._c)
            new._root = root
            new._len = self._len - 1
            return new
        else:
            return self

    def get(self, x):
        self._check_x(x)
        found, data = PVeb._get(self._root, x - self._lb, self._u, self._threshold)
        if found:
            return data
        else:
            raise KeyError(f"Key: {x} not found!")

    def contains(self, x):
        self._check_x(x)
        return PVeb._get(self._root, x - self._lb, self._u, self._threshold)[0]

    def successor(self, x):
        self._check_x(x)
        ret = PVeb._successor(self._root, x - self._lb, self._u, self._threshold)
        if ret != None:
            k, v = ret
            return (k + self._lb), v

    def predecessor(self, x):
        self._check_x(x)
        ret = PVeb._predecessor(self._root, x - self._lb, self._u, self._threshold)
        if ret != None:
            k, v = ret
            return (k + self._lb), v

    def min(self):
        return (
            ((self._lb + self._root.min), self._root.min_data) if self._root else None
        )

    def extractMin(self):
        if self._root:
            root = PVeb._pop(self._root, self._root.min, self._u, self._threshold)
            new = PVeb(lb=self._lb, ub=self._ub, c=self._c)
            new._root = root
            new._len = self._len - 1
            return new
        else:
            return self

    def max(self):
        return (
            ((self._lb + self._root.max), self._root.max_data) if self._root else None
        )

    def extractMax(self):
        if self._root:
            root = PVeb._pop(self._root, self._root.max, self._u, self._threshold)
            new = PVeb(lb=self._lb, ub=self._ub, c=self._c)
            new._root = root
            new._len = self._len - 1
            return new
        else:
            return self

    def iter_items(self, reverse=False):
        if self._root:
            f = self.predecessor if reverse else self.successor
            kv = self.max() if reverse else self.min()

            while kv:
                yield kv
                kv = f(kv[0])

    def iter_keys(self, reverse=False):
        return (k for k, _ in self.iter_items(reverse=reverse))

    def items(self, reverse=False):
        return list(self.iter_items(reverse=reverse))

    def keys(self, reverse=False):
        return list(self.iter_keys(reverse=reverse))


Mapping.register(PVeb)
Hashable.register(PVeb)
