### Configurations ###
######################

# Computation Settings
solutions_limit = -1 # how many solution you wish to find
interactions_limit = 17 # limit of optional interactions.
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
          "find_minimal_combination": ( "SimplestModel.txt", # establised model
                                        "UltimateConstrains.txt" ),
          "find_minimal_model":
              ( "PearsonThreshold792WithOct4Sox2Interaction.txt",
                "UltimateConstrains.txt" )}
MODEL, EXP = INPUT['ABCD_test']

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

    if __debug__:
        print '>> Species codes:'; pprint(code)
        print '>> Defined Iteractions:'; pprint(defI)
        print '>> Optional Interactions:'; pprint(optI)

    bitlen = len(species)
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
        f_[s] = [makeFunction(acts, reps, l, A_[s], R_[s]) for l in logics[s]]

    solver = Solver()
    # setup functions in solver
    bv = BitVec('bv', len(species)) # used in ForAll expression. Used locally
    F = Function('F', bv.sort(), bv.sort()) # just declaration
    for s in species:
        c = code[s]
        for l in range(L_[s].size()):
            solver.add(Implies(
                Extract(l,l,L_[s]) == 1, # when l is selected
                ForAll(bv, (Extract(c,c,F(bv))==1) == f_[s][l](bv))))
        if __debug__: print '>> Set function for %s: %s'%(s,solver.check())
        
    # define Transition ralationship. it is a macro.
    T = lambda qo, qn: qn == F(qo)
    
    # Apply modeling constrains
    for s in species:
        c = code[s]
        # defined activators and repressors must be selected
        defact, defrep = defI[c]
        if defact: solver.add(1 == Extract(len(defact)-1, 0, A_[s]))
        if defrep: solver.add(1 == Extract(len(defrep)-1, 0, R_[s]))
        # only one logic is selected
        logic_i = range(L_[s].size())
        solver.add(1 == ~(Any([ Extract(i,i,L_[s]) & Extract(j,j,L_[s]) \
                                for i in logic_i for j in logic_i if i != j])))
        # must select one logic
        solver.add(L_[s] != 0)
        
        if __debug__:
            print '>> Constraints %s:\tAct = %s,\tRep = %s,\tLog = %s.\t(%s)' \
                %(s, A_[s] and A_[s].sort() or None,
                  R_[s] and R_[s].sort() or None,
                  L_[s].sort(), solver.check())
        
    # apply experimental constrains
    for exp in exps:
        if __debug__: print '>> Adding %s ... '%exp
        # build path
        path = [BitVec(exp + '_%d'%t, bitlen) for t in range(STEP)]
        solver.add(*[ T(path[t], path[t+1]) for t in range(STEP-1) ])
        # add constrains
        for t, conditions in exps[exp]:
            for cond in conditions:
                for s, value in states[cond]:
                    c = code[s]
                    solver.add( Extract(c,c,path[t]) == value )
    print ">> Constrains established."

    count = 0
    allAR = list(A_.values()) + list(R_.values())
    while solver.check() == sat:
        count += 1
        # make sure all A_[s] and R_[s] are specified
        m = solver.model()
        solver.push()
        for b in allAR:
            if b and not m[b]:
                solver.add( b == 0)
        assert solver.check() == sat
        m = solver.model() # update the solution
        # print out
        print ">> Solution %d: "%count
        printModel(m, A_, R_, L_, species, code, inters)
        if count == solutions_limit: break
        # find different solutions (with different selections of interactions)
        # at least one species have distinct interactions
        solver.pop()
        solver.add(Or([ b != m[b] for b in allAR if b ]))
                      
    if count == 0: print 'No solution found.'

if __name__ == '__main__':
    main()
