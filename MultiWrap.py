#!/usr/bin/env python

import mapfn
import resultfn

# def FunWrap(thisfun):
#     def FunOut(args):
#         return thisfun(*args)
#     return FunOut

def MakeWrap(thisfun,constarg):
    def NewFun(iterarg):
        return thisfun(*(constarg+(iterarg,)))
    return NewFun

def MakeWrapD2(thisfun,constarg):
    def NewFun(iterarg):
        return thisfun(*(constarg+(iterarg[0],iterarg[1])))
    return NewFun


def makeContextFunctions(*fns):
    for f in fns:
        newName = f.__name__
        def mapper(x):
            return f(*x)
        mapper.__name__ = newName
        mapper.__module__ = 'mapfn'
        setattr(mapfn, newName, mapper)
        f.mapper = mapper

        def resultFn(x):
            print "%s returned result %s" % (f.__name__, x)
        resultFn.__name__ = newName
        resultFn.__module__ = resultfn
        setattr(resultfn, newName, resultFn)
        f.resultFn = resultFn


# def MakeCWrap(thisfun,constargs):
#     funout = MakeWrap(thisfun,constargs)
#     makeContextFunctions(funout)
#     return funout

# def MakeCWrapD2(thisfun,constargs):
#     funout = MakeWrapD2(thisfun,constargs)
#     makeContextFunctions(funout)
#     return funout
