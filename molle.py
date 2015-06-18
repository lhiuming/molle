### Configurations ###
######################

# Computation Settings
SOLUTIONS_LIMIT = 0 # how many solution you wish to find
INTERACTIONS_LIMIT = 17 # limit of optional interactions. 17 for minimal model

# Input and output files
PREFIX = "examplefiles/"
INPUT = { 'ABCD_test': ( "SimpleFourComponentModel.txt",
                         "CertainInteractionRequired.txt" ), # not true
          'ABCD_kofe': ( "four_modified.txt",
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
MODEL, EXP = INPUT['find_min_inter']

# Model Configurations
STEP = 20 # trajactory length

### Modelling #########################################
#######################################################

from z3 import *
from utility import *
from pprint import pprint
from time import time, localtime, strftime

def main(solver, solutions_limit=10, interactions_limit=0,
         mfile=MODEL, efile=EXP, debug=False, detail=True, output = True):
    startt = time()
    
    # reading files
    modelFile = open(PREFIX + mfile, 'r')
    (species, code, logics, kofe, defI, optI) = readModel(modelFile)
    modelFile.close()
    expFile = open(PREFIX + efile, 'r')
    (exps, states) = readExp(expFile)
    expFile.close();

    # convinient variables
    bitlen = len(species)
    kos, fes = kofe['KO'], kofe['FE']
    
    if debug:
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
        kofe_index = (s in kos and kos.index(s)+1,
                      s in fes and fes.index(s)+1)
        f_[s] = [makeFunction(acts, reps, kofe_index, l, A_[s], R_[s])\
                 for l in logics[s]]


    # 1. Modeling constrains
    if detail: print '>> #1 Adding modeling constrains: '
    for s in species:
        c = code[s]
        # INTER: defined activators and repressors must be selected
        defact, defrep = defI[c]
        solver.add([1 == Extract(i,i,A_[s]) for i in range(len(defact))])
        solver.add([1 == Extract(i,i,R_[s]) for i in range(len(defrep))])
        # LOGIC: only one logic is selected
        logic_i = range(L_[s].size())
        solver.add(0 == Any([ Extract(i,i,L_[s]) & Extract(j,j,L_[s])
                              for i in logic_i for j in logic_i if i != j]))
        # LOGIC: must select one logic
        solver.add(L_[s] != 0)
        if debug:
            assert solver.check() == sat
            print '>> \tConstraints %s:\tAct=%s,\tRep=%s,\tLog=%s(%s).' \
                %(s, A_[s] and A_[s].sort() or None,
                  R_[s] and R_[s].sort() or None,
                  L_[s].sort(), isExpOf2(solver.model()[L_[s]]) or
                  bin(solver.model()[L_[s]].as_long()))
            
    # 2. Interactions limit (only for optional interactions)
    opts = [] # all optional interactions
    if interactions_limit:
        for s in species:
            c = code[s]; actn, repn = map(len,optI[c]) # nums of ats and reps
            # fill the allOpt list with all optional interactions
            if A_[s]:
                opts.extend([Extract(i,i,A_[s])
                             for i in range(A_[s].size()-actn, A_[s].size())])
            if R_[s]:
                opts.extend([Extract(i,i,R_[s]) \
                             for i in range(R_[s].size()-repn, R_[s].size())])
        # constraints that all selected nums of inters are less than limit
        solver.add(ULE(sum([ZeroExt(6, b) for b in opts]),
                       interactions_limit))
        if detail:
            print '>> #2 Interactions limit %d ADDED. %s' \
                %(interactions_limit,solver.check())
    else:
        if detail: print '>> #2 Interactions limit: NOT SET.'
        
    
    # 3. Setup functions in solver
    #bv = BitVec('bv', bitlen) # used in ForAll expression. Used locally
    #ko = BitVec('ko', len(kos) or 1)
    #fe = BitVec('fe', len(fes) or 1)
    #F = Function('F', bv.sort(), ko.sort(), fe.sort(), bv.sort()) # declare
    #for s in species:
    #    c = code[s]
    #    for l in range(L_[s].size()):
    #        solver.add(Implies(Extract(l,l,L_[s]) == 1, # when l is selected
    #                           ForAll([bv, ko, fe],
    #                                  (Extract(c,c,F(bv, ko, fe))==1) == \
    #                                  f_[s][l](bv, ko, fe))))
    #    if __debug__: print '>> Set function for %s: %s'%(s, solver.check())
    if detail: print '>> #3 Setup forall style function-definition: SKIP.'


    # 4. Define Transition ralationship. It is like a macro.
    T = lambda q_old, q_new, ko, fe: \
        And([ And([Implies(Extract(l,l, L_[s]) == 1, # if logic l is selected
                       (Extract(code[s],code[s],q_new)==1) == \
                        f_[s][l](q_old, ko, fe)) # update value
                  for l in range(L_[s].size())])
               for s in species ])
    bunchT = lambda qs, ko, fe: \
             And([ And([ Implies(Extract(l,l,L_[s]) == 1,
                               And([ (Extract(code[s],code[s],qs[t+1])==1) == \
                                     f_[s][l](qs[t], ko, fe)
                                     for t in range(STEP-1) ]))
                         for l in range(L_[s].size())])
                   for s in species ])
    if detail: print '>> #4 Transition relation macro: SET.'

    # 5. Experimental constrains
    if detail: print '>> #5 Applying experimental constraints:'
    total_exp = len(exps)
    count_exp = 0
    for exp in exps:
        count_exp += 1
        # build path
        KO = BitVec(exp + '_KO', len(kos) or 1)
        FE = BitVec(exp + '_FE', len(fes) or 1)
        path = [BitVec(exp + '_%d'%t, bitlen) for t in range(STEP)]
        #solver.add(*[ T(path[t], path[t+1], KO, FE) for t in range(STEP-1) ])
        solver.add(bunchT(path, KO, FE)) # a little bit faster in solving
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
                        c = code[s]
                        solver.add( Extract(c,c,path[t]) == value )
        if detail:
            print '>> \t %02d/%d %s added...'%(count_exp, total_exp, exp)
    if detail: print ">> All Constrains established."

    # Now get the solutions !
    solvingt = time() # just for timing
    print '>> ' + '-'* 76
    print '>> Start solving: %s'%strftime("%d %b %H:%M", localtime(startt))
    print '>> '
    
    count = 0
    allAR = [b for b in list(A_.values()) + list(R_.values()) if b]
    while solver.check() == sat:
        count += 1
        # make sure all A_[s] and R_[s] are specified
        solver.push() # push will somehow change the solution
        solver.check(); m = solver.model()
        solver.add([ b == 0 for b in allAR if not m[b] ]) # must hola a value
        solver.add([ b == m[b] for b in allAR if m[b] ]) # keep the  value
        solver.check(); m = solver.model() # update the solution
        solver.pop()
        # print out
        if output:
            print ">> Solution %d: "%count
            printModel(m, A_, R_, L_, species, code, inters, logics,
                       config = True, model = True)
        if count == solutions_limit: break
        # find different solutions (with different selections of interactions)
        # at least one species have distinct interactions
        solver.add(Or([ b != m[b] for b in allAR]))
        
    endt = time()
    print '>> %d solutions found. \n>> '%count
    print '>> End solving: %s.'%strftime("%d %b %H:%M",localtime(endt))
    print '>> ' + '-'* 76
    if endt - startt > 600:
        solving = '%.1f min'%((endt - solvingt) / 60)
        total = '%.1f min'%((endt - startt) / 60)
    else:
        solving = '%.2f s'%(endt - solvingt)
        total = '%.2f s'%(endt - startt)
    print '>> Solving duration:\t%s'%solving
    print '>> Total duration:\t%s'%total
    print '>> ' + '-' * 3 + ' Finished. ' + '-' * 3

if __name__ == '__main__':
    s = Solver()
    main(s, SOLUTIONS_LIMIT, INTERACTIONS_LIMIT, MODEL, EXP)
