# temporary utility file

from itertools import combinations
from z3 import *
from utility import _create_rule
from operator import and_

def sortedInters(inters, code):
    ''' sorted the inters into a dict, where keyvalue is a couple of number
    tuples.'''
    d = dict([(c, ([], [])) for c in code.values()])
    for i in inters:
        f, t = i[:2]
        idx = (0, 1)[ i[2]=='negative' ]
        d.setdefault(code[t], ([], []))[idx].append(code[f])
    return d

def _sublist_generator(l, high=None, low=0):
    if not high: high = len(l)
    for r in range(low, high + 1):
        g = combinations(l, r)
        for sl in g: yield sl

def getInterCombi(defs, opts):
    ''' Return a list of all iteraction combinations, represented by a couple
    of tuples, containing activators and represors.
    '''
    defAct, defRep = defs
    optAct, optRep = opts
    return [(tuple(defAct + list(act)), tuple(defRep + list(rep))) \
            for act in _sublist_generator(optAct) \
            for rep in _sublist_generator(optRep)]

def _concat(l):
    if not l:return BitVecVal(0, 1)
    if len(l) == 1: return l[0]
    return Concat(l)

def _create_rule(num, act, rep):
    if num == -1: return False
    if num < 2 and rep: return False
    if num > 15 and act: return False
    
    actB = _concat(act)
    repB = _concat(rep)

    if num==0: return ~actB == 0
    elif num==1: return actB > 0
    elif num<4: return And(actB > 0, repB == 0)
    elif num<6: return And(~actB == 0, ~repB > 0)
    elif num<8: return And(actB > 0, ~repB > 0)
    elif num<10: return ~actB == 0
    elif num<12: return Or(~actB == 0, And(actB > 0, repB == 0))
    elif num<14: return Or(~actB == 0, And(actB > 0, ~repB > 0))
    elif num<16: return actB > 0
    elif num==16: return And(repB > 0, ~repB > 0)
    elif num==17: return repB == 0
    
def makeFunction(inter, logic_num):
    ''' make a function that takes a q, and return a coresponding z3 expr.'''
    return lambda q: _create_rule(logic_num,
                                  [Extract(i, i, q) for i in inter[0]],
                                  [Extract(i, i, q) for i in inter[1]])

if __name__ == "__main__":
  if __debug__:
    print ">> testing sortedInters():"
    print sortedInters([['a', 'b', 'positive'], ['b', 'c', 'negative']],
                       {'a':0, 'b':1, 'c':2})

    print ">> testing getInterCombi():"
    testDef = (['a1', 'a2'], ['r1', 'r2'])
    testOpt = (['aaa'], [])
    print getInterCombi(testDef, testOpt)

    print ">> testing makeFunction():"
    testq = BitVec('testq', 8)
    testI = ( (1, 3), (6, 7) )
    print makeFunction(testI, 1)(testq)
    print makeFunction(testI, 16)(testq)
    print makeFunction(testI, 7)(testq)
