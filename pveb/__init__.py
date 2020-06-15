from ._pveb import PVeb

_PVEB_CACHE = {}


def pveb(lb=0, ub=(1 << 32) - 1, c=100):
    """This module returns a persistent van-emde-boas tree
    
    Keyword Arguments:
        lb {int} -- [lb of the van-emde-boas tree] (default: {0})
        ub {int} -- [ub of the van-emde-boas tree] (default: {(1 << 32)-1})
        c {int} -- [
            parameter to define van-emde-boas height threshold, higher the c,
            lower the height, hence slower can be the search. Choose wisely
        ] (default: {100})
    """
    z = lb, ub, c
    if z in _PVEB_CACHE:
        return _PVEB_CACHE[z]
    else:
        _PVEB_CACHE[z] = v = PVeb(lb=lb, ub=ub, c=c)
        return v
