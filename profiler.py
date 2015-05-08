import cProfile
import molle

cProfile.run('molle.main()', 'mollestats')

import pstats
p = pstats.Stats('mollestats')
p.strip_dirs().sort_stats('time').print_stats()
