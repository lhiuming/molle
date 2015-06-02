### Configurations ###
######################

# Running Settings
debug = 1 # debugging level

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

    # encoding all devices/modules/functions
    I_ = [] # of Interaction-selecting BitVec for species
    L_ = [] # of Logic-selecting BitVec for species
    iNum = [] # total iteraction-combination number for species
    f_ = [] # devices/modules/functions
    for s in species:
        c = code[s]
        # make function for all combinations
        f_[c] = []
        i = 0
        for i, inter in enumerate(generateInterCombi(defI[c], optI[c])):
            f_[c][i] = dict([ (l, makeFunction(inter, l)) for l in logics[s] ])
        iNum[c] = i + 1
        # create interaction-selecting BitVec
        if i > 0:
            I_[c] = BitVec('I_' + s, iNum[c])
        else: # no optional interaction
            I_[c] = BitVecVal(0, 1)
        # create logic-selecting BitVec
        L_[c] = BitVec('L_' + s, len(logcis[s]))

    # define Transition ralationship
    T = lambda qo, qn:\
        And(*[ Extract(c, c, qn) == (
                   Extract(i, i, I_[c]) &
                   Extract(l, l, L_[c]) &
                   f_[c][i][l](q0))
              for c in code.values()])
            
