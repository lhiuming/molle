from molle import ABN
from utility import *
from time import time, localtime, strftime
import sys

PREFIX = "examplefiles/"
INPUT = { 'ABCD_test': ( "SimpleFourComponentModel.txt",
                         "CertainInteractionRequired.txt" ), # not true
          'ABCD_kofe': ( "four_modified.txt",
                         "four_constrainst.txt" ),
          'logics_range_test': ( "four_logic_modified.txt",
                                "four_constrainst.txt" ),   
          'ABCD_nosolution': ("SimpleFourComponentModel.txt",
                              "NoSolutionsPossible.txt" ),
          'minimal_test': ( "custom.txt", # established model%combination
                            "UltimateConstrains.txt" ),
          "find_min_logic": ( "simplestmodel.txt", # establised inters
                              "UltimateConstrains.txt" ),
          'find_min_inter': ( "simplestlogic.txt", # established logics
                              "UltimateConstrains.txt" ),
          "find_minimal_model":
              ( "PearsonThreshold792WithOct4Sox2Interaction.txt",
                "UltimateConstrains.txt" )}

if __name__ == '__main__':
    # read commond line argvs
    nonpar = filter(lambda x: '-' != x[0], sys.argv[1:3])
    if len(nonpar) == 2:
        mpath, epath = nonpar
        problem = mpath + ', ' + epath
    elif len(nonpar) == 1:
        problem = nonpar[0]
        mpath, epath = [PREFIX + name for name in INPUT[problem]]
    else:
        print "No enough arguments. ending."; quit()
    output = '-o' in sys.argv
    debug = '-d' in sys.argv
    verbose = '-v' in sys.argv or debug
    mail = '-m' in sys.argv
    if mail:
        i = sys.argv.index('-m') + 1
        addr, pw = sys.argv[i:i+2]

    # start timing
    startt = time()
    print '>> '+  "Abstract Network Modeling: molle " + '-' * 48
    print '>> Modeling for %s'%problem
    print '>> Start program: %s'%strftime("%d %b %H:%M", localtime(startt))    

    # reading inputs
    modelFile = open(mpath, 'r')
    expFile = open(epath, 'r')
    model = ABN(modelFile, expFile)
    modelFile.close(); expFile.close()

    # build the model
    model.build(detail=verbose, debug=debug)
    print ">> All Constrains established. (takes %s)" %conv_time(time()-startt)

    # get solutions
    print '>> ' + '- '* 15 # seperator
    solvingt = lastt = time() # just for timing
    print '>> Start solving: %s'%strftime("%d %b %H:%M",localtime(solvingt))
    count = 0
    for solution in model.solve():
        count += 1
        print ">> Solution %d: (takes %.1f min)"%(count,(time() - lastt)/60)
        lastt = time() # update time for lasted model
        if output: solution.output()
    endt = time()

    # concluding
    print '>> %d solutions found. (takes %s to end)' \
        %(count, conv_time(endt - lastt))
    print '>> End solving: %s.'%strftime("%d %b %H:%M",localtime(endt))
    print '>> ' + '- ' * 15
    print '>> Solving duration:\t%s'%conv_time(endt - solvingt)
    total = conv_time(endt - startt)
    print '>> Total duration:\t%s'%total
    if mail:
        mailMe(addr, pw,
               'Solutions number:\t%d\nTotal duration:\t%s'%(count, total),
               "Computation Fnishied for '%s'."%problem)
    print '>> ' + '-' * 9 + ' Finished. ' + '-' * 9





