from z3 import Solver
import molle

s = Solver()
m,e = molle.INPUT['ABCD_kofe']
print '** 1 ** TESTING ABCD_kofe: '
molle.main(s, 0, 0, m, e, debug = False, detail = False, output = False)

s.reset()
m,e = molle.INPUT['minimal_test']
print '** 2 ** TESTING minimal_test: '
molle.main(s, 0, 17, m, e, debug = False, detail = False, output = False)
