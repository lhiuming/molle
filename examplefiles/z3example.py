# "Mapping Functions on Arrays" Section on tutorial
#
# Notes: to use identifier _, which is not included in the python API, we use
#        parse_smt2_string to assert expressions like (_ map and).

from z3 import *

Set = ArraySort(IntSort(), BoolSort()) # set sort
a, b, c = Consts('a b c', Set) # three sets
x = Int('x')
s = Solver()

s.push()
s.add(parse_smt2_string('''
(assert (not (= ((_ map and) a b)
                ((_ map not) ((_ map or) ((_ map not) b)
                                         ((_ map not) a))))))''',
                        decls={'a':a, 'b':b}))
print s.check()

s.pop()
s.push()
s.add(parse_smt2_string('''
(assert (and (select ((_ map and) a b) x)
             (not (select a x))))''', decls={'a':a, 'b':b, 'x':x}))
print s.check()

s.pop()
s.add(parse_smt2_string('''
(assert (and (select ((_ map or) a b) x)
             (not (select a x))))''', decls={'a':a, 'b':b, 'x':x}))
print s.check() # sat
print s.model().sexpr()
s.add(And(Not(Select(b, x))))
print s.check()


