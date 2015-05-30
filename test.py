from z3 import *
from utility import readModel,  readExp
from utilitytmp import *

PREFIX = "examplefiles/"
INPUT = { 'ABCD_test': ( "SimpleFourComponentModel.txt",
                         "CertainInteractionRequired.txt" ),}
MODEL, EXP = INPUT['ABCD_test']
modelFile = open(PREFIX + MODEL, 'r')
(logics, kofe, defInters, optInters) = readModel(modelFile); modelFile.close()
expFile = open(PREFIX + EXP, 'r')
(exps, states) = readExp(expFile); expFile.close();

species = comps = logics.keys()
code = dict([(c, n) for (n, c) in enumerate(species)])
print code
print logics

# D = BitVec('D', diveceNum) # encoding the set of selected devices
# Q = BitVecSort(compNum) # state sort

## encoding all devices/functions
I_ = dict();L_ = dict();interNum = dict();f_ = {}
defIs = sortedInters(defInters, code)
optIs = sortedInters(optInters, code)
for s in species:
    c = code[s] # a unique integral code for the specie
    inter_combis = getInterCombi(defIs[c], optIs[c])
    interNum[s] = len(inter_combis)
    f_[s] = dict() # hold concrete functions
    for i, inter in enumerate(inter_combis):
        f_[s][i] = dict([(l, makeFunction(inter, l)) for l in logics[s]])
    if interNum[s]:
        I_[s] = BitVec('I_' + s, interNum[s]) # BitVec for selecting inters
    else: I_[s] = BitVecVal(0, 1)
    L_[s] = BitVec('L_' + s, len(logics[s])) # BitVec for selecting logic

## define transition
T = lambda qo, qn: And(*[ (1 == Extract(code[s], code[s], qn)) == \
                            Or(*[ And(1 == (Extract(i, i, I_[s]) & \
                                            Extract(l, l, L_[s])),
                                      f_[s][i][l](qo)) \
                                  for i in range(interNum[s]) \
                                  for l in logics[s]] ) \
                            for s in species])

## applyr constraints
solver = Solver()
# only one interaction-combination and one logic for a specie
for s in species:
    if interNum[s] > 1:
        solver.add(1 == reduce(and_,
                          [~(Extract(i, i, I_[s]) & Extract(j, j, I_[s])) \
                           for i in range(interNum[s]) \
                           for j in range(interNum[s]) if i != j]))
    if len(logics[s]) > 1:
        solver.add(1 == reduce(and_,
                          [~(Extract(i, i, L_[s]) & Extract(j, j, L_[s])) \
                           for i in logics[s] for j in logics[s] if i != j]))

#print solver.check()
#if solver.check() == sat:
#    print solver.model().sexpr()

#print exps
#print states

bitlen = len(species)
for exp in exps:
    print exp
    # build a trajectory/path of 20 steps
    path = [BitVec(exp + '_%d'%t, bitlen) for t in range(20)]
    solver.add(And(*[ T(path[t], path[t+1]) for t in range(19) ]))

if solver.check() == sat:
    print solver.model()
