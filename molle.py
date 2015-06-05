### Configurations ###
######################

# Computation Settings
solutions_limit = 0 # how many solution you wish to find
interactions_limit = 16 # limit of optional interactions.
                        # set to 17 for the minimal pluripotency model

# Input and output files
OUTPUT = "solutions.out" # file name for solution output
PREFIX = "examplefiles/"
INPUT = { 'ABCD_test': ( "SimpleFourComponentModel.txt",
                         "CertainInteractionRequired.txt" ), # not true
          'ABCD_nosolution': ("SimpleFourComponentModel.txt",
                              "NoSolutionsPossible.txt" ),
          'minimal_test': ( "custom.txt", # established model%combination
                            "UltimateConstrains.txt" ),
          'minimal_test2': ( "custom2.txt", # established model%combination
                            "UltimateConstrains.txt" ),
          "find_minimal_combination": ( "SimplestModel.txt", # establised model
                                        "UltimateConstrains.txt" ),
          "find_minimal_model":
              ( "PearsonThreshold792WithOct4Sox2Interaction.txt",
                "UltimateConstrains.txt" )}
MODEL, EXP = INPUT['minimal_test2']

# Model Configurations
STEP = 20 # trajactory length


### Modelling #########################################
#######################################################

from z3 import *
from utility import *
from pprint import pprint

def main():
    # reading files
    modelFile = open(PREFIX + MODEL, 'r')
    (species, code, logics, kofe, defI, optI) = readModel(modelFile)
    modelFile.close()
    expFile = open(PREFIX + EXP, 'r')
    (exps, states) = readExp(expFile)
    expFile.close();

    bitlen = len(species)
    kos, fes = kofe['KO'], kofe['FE']
    
    if __debug__:
        print '>> Species codes:'; pprint(code)
        print '>> Defined Iteractions:'; pprint(defI)
        print '>> Optional Interactions:'; pprint(optI)
        print '>> KO: ', kos
        print '>> FE: ', fes

    # encoding all devices/modules/functions
    A_ = {} # of activator/activating-interaction selection BitVec
    R_ = {} # of repressor/repressing-interaction selection BitVec
    L_ = {} # of Logic-selecting BitVec for speciew
    f_ = {} # devices/modules/functions
    inters = {} # record all avalible acts and reps for every specie
    for s in species:
        c = code[s] # numeral index
        acts = optI[c][0] + defI[c][0] # Concat is from left to right
        reps = optI[c][1] + defI[c][1]
        inters[c] = (acts, reps)
        # creating Acr and Rep selecting BitVec
        if acts: A_[s] = BitVec('Act_' + s, len(acts))
        else: A_[s] = None
        if reps: R_[s] = BitVec('Rep_' + s, len(reps))
        else: R_[s] = None
        # create logic-selecting BitVec
        L_[s] = BitVec('Logic_' + s, len(logics[s]))
        # make the functions
        kofe_index = (s in kos and kos.index(s)+1, s in fes and fes.index(s)+1)
        f_[s] = [makeFunction(acts, reps, kofe_index, l, A_[s], R_[s])\
                 for l in logics[s]]

    solver = Solver()
    # setup functions in solver
    bv = BitVec('bv', bitlen) # used in ForAll expression. Used locally
    ko = BitVec('ko', len(kos) or 1)
    fe = BitVec('fe', len(fes) or 1)
    F = Function('F', bv.sort(), ko.sort(), fe.sort(), bv.sort()) # declaration
    for s in species:
        c = code[s]
        for l in range(L_[s].size()):
            solver.add(Implies(
                Extract(l,l,L_[s]) == 1, # when l is selected
                ForAll([bv, ko, fe], (Extract(c,c,F(bv, ko, fe))==1) == \
                       f_[s][l](bv, ko, fe))))
        if __debug__: print '>> Set function for %s: %s'%(s, solver.check())

    # Apply modeling constrains
    for s in species:
        c = code[s]
        # defined activators and repressors must be selected
        defact, defrep = defI[c]
        solver.add([1 == Extract(i,i,A_[s]) for i in range(len(defact))])
        solver.add([1 == Extract(i,i,R_[s]) for i in range(len(defrep))])
        # only one logic is selected
        logic_i = range(L_[s].size())
        solver.add(1 == ~(Any([ Extract(i,i,L_[s]) & Extract(j,j,L_[s]) \
                                for i in logic_i for j in logic_i if i != j])))
        # must select one logic
        solver.add(L_[s] != 0 )
        if __debug__:
            print '>> Constraints %s:\tAct = %s,\tRep = %s,\tLog = %s. (%s)' \
                %(s, A_[s] and A_[s].sort() or None,
                  R_[s] and R_[s].sort() or None,
                  L_[s].sort(), 'nocheck' or solver.check())

            # define Transition ralationship. it is a macro.
    T = lambda q_old, q_new, ko, fe: q_new == F(q_old, ko, fe)

    # interactions limit
    allOpt = []
    for s in species:
        c = code[s]; actn, repn = map(len,optI[c])
        if A_[s]:
            allOpt.extend([Extract(i,i,A_[s]) \
                       for i in range(A_[s].size()-actn, A_[s].size())])
        if R_[s]:
            allOpt.extend([Extract(i,i,R_[s]) \
                       for i in range(R_[s].size()-repn, R_[s].size())])
    solver.add(ULE(sum([ZeroExt(6, b) for b in allOpt]), interactions_limit))
    print '>> Interactions limit added. %s'%solver.check()
        
    # apply experimental constrains
    for exp in exps:
        if __debug__: print '>> Adding %s ... '%exp
        # build path
        KO = BitVec(exp + '_KO', len(kos) or 1)
        FE = BitVec(exp + '_FE', len(fes) or 1)
        path = [BitVec(exp + '_%d'%t, bitlen) for t in range(STEP)]
        solver.add(*[ T(path[t], path[t+1], KO, FE) for t in range(STEP-1) ])
        # add constrains
        for t, conditions in exps[exp]:
            for cond in conditions:
                for s, value in states[cond]:
                    sl = s.split('_')
                    s = sl[-1]
                    if sl[0] == 'KO':
                        c = kos.index(s)
                        solver.add( Extract(c,c,KO) == value)
                    elif sl[0] == 'FE':
                        c = fes.index(s)
                        solver.add( Extract(c,c,FE) == value)
                    else:
                        c = code[sl[-1]]
                        solver.add( Extract(c,c,path[t]) == value )
    print ">> Constrains established."

    count = 0
    allAR = [b for b in list(A_.values()) + list(R_.values()) if b]
    while solver.check() == sat:
        count += 1
        # make sure all A_[s] and R_[s] are specified
        solver.push() # push will somehow change the solution
        solver.check(); m = solver.model()
        solver.add([ b == 0 for b in allAR if not m[b] ])
        solver.check(); m = solver.model() # update the solution
        solver.pop()
        # print out
        print ">> Solution %d: "%count
        printModel(m, A_, R_, L_, species, code, inters,\
                   config = True, model = True)
        if count == 3: print m
        if count == solutions_limit: break
        # find different solutions (with different selections of interactions)
        # at least one species have distinct interactions
        solver.add(Or([ b != m[b] for b in allAR]))                      
    if count == 0: print 'No solution found.'

if __name__ == '__main__':
    main()
