### Configurations ###
######################

# Computation Settings
solutions_limit = 0 # how many solution you wish to find
interactions_limit = 17 # limit of optional interactions.
                        # set to 17 for the minimal pluripotency model

# Input and output files
OUTPUT = "solutions.out" # file name for solution output
PREFIX = "examplefiles/"
INPUT = { 'ABCD_test': ( "SimpleFourComponentModel.txt",
                         "CertainInteractionRequired.txt" ),
          'ABCD_nosolution': ("SimpleFourComponentModel.txt",
                              "NoSolutionsPossible.txt" ),
          'minimal_test': ( "custom.txt", # established model%combination
                            "UltimateConstrains.txt" ),
          "find_minimal_combination": ( "SimplestModel.txt", # establised model
                                        "UltimateConstrains.txt" ),
          "find_minimal_model":
              ( "PearsonThreshold792WithOct4Sox2Interaction.txt",
                "UltimateConstrains.txt" )
        }
MODEL, EXP = INPUT['ABCD_test']

# Model Configurations
STEP = 20 # trajactory length
use_compact = False # use compact list of allowed function


### Modelling #########################################
#######################################################

from z3 import *
from utility import *
from pprint import pprint

#from multiprocessing import Process, Pipe
#from Queue import Queue

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
    # I_ = {} # of Interaction-selecting BitVec for species
    A_ = {} # activator selecting
    R_ = {} # repressor selecting
    L_ = {} # of Logic-selecting BitVec for species
    f_ = {} # devices/modules/functions
    # make function for all logics
    for s in species:
        c = code[s]
        acts = defI[c][0] + optI[c][0] # defined elements is before optionals
        reps = defI[c][1] + optI[c][1]
        # creating A/R selecting BitVec
        if acts: A_[s] = BitVec('Act_' + s, len(acts))
        else: A_[s] = BoolVal(False)
        if reps: R_[s] = BitVec('Rep_' + s, len(reps))
        else: R_[s] = BoolVal(False)
        # create logic-selecting BitVec
        L_[s] = BitVec('Logic_' + s, len(logics[s]))
        # make the functions
        f_[s] = [makeFunction(acts, reps, l, A_[s], R_[s]) for l in logics[s]]

    solver = Solver()

    # setup functions in solver
    bv = BitVec('bv', len(species)) # used in ForAll expression
    F = Function('F', bv.sort(), bv.sort()) # just declaration
    for s in species:
        c = code[s]
        for l in range(L_[s].size()):
            solver.add(Implies(
                Extract(l,l,L_[s]) == 1,
                ForAll(bv, (Extract(c,c,F(bv))==1) == f_[s][l](bv))))
        if __debug__: print '>> Set function for %s: %s'%(s,solver.check())
        
    # define Transition ralationship. it is a macro.
    T = lambda qo, qn: qn == F(qo)
    
    ## apply model constrains
    # only one interaction-combination and one logic for a gene/species
    for s in species:
        c = code[s]
        # defined activators and repressors must be valid
        defact = defI[c][0]; defrep = defI[c][1]
        if defact: solver.add(1 == Extract(len(defact)-1, 0, A_[s]))
        if defrep: solver.add(1 == Extract(len(defrep)-1, 0, R_[s]))
        # only one logic is selected
        solver.add(1 == ~(any([ Extract(i,i,L_[s]) & Extract(j,j,L_[s]) \
                                for i in range(L_[s].size()) \
                                for j in range(L_[s].size()) if i != j])))
        # must select one logic
        solver.add(L_[s] != 0)
        
        if __debug__: print '>> Constrainting %s: A = %s, R = %s, L = %s, %s' \
           %(s, A_[s].sort(), R_[s].sort(), L_[s].sort(), solver.check())
        
    # apply experimental constrains
    for exp in exps:
        # build path
        if __debug__: print '>> Adding %s... '%exp
        path = [BitVec(exp + '_%d'%t, bitlen) for t in range(STEP)]
        solver.add(And(*[ T(path[t], path[t+1]) for t in range(STEP-1) ]))
        # add constrains
        for t, conditions in exps[exp]:
            for cond in conditions:
                for species, value in states[cond]:
                    c = code[species]
                    solver.add( Extract(c,c,path[t]) == value )
    if solver.check() == sat:
        m = solver.model()
        print m
    else: print 'No solution found.'

if __name__ == '__main__':
    main()
