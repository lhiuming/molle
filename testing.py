from z3 import Solver
import molle

s = Solver()
print '** 1 ** TESTING ABCD_kofe: '
molle.main(s, 'ABCD_kofe', 0, 0, debug = False, detail = False, output = False)

s.reset()
print '** 2 ** TESTING minimal_test: '
molle.main(s, 'minimal_test', 0, 17, debug = False, detail = False, output = False)
