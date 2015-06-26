### Configurations ###
######################

# Computation Settings
STEP = 20 # trajactory length

# Input and output files
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

### Modelling #########################################
#######################################################

from z3 import *
from utility import *
from pprint import pprint
from time import time, localtime, strftime

def main(solver, problem, solutions_limit=10, interactions_limit=0,
         debug=False, detail=True, output = True):
    startt = time()
    print '>> Start program: %s'%strftime("%d %b %H:%M", localtime(startt))    

    # reading files
    mfile, efile = INPUT[problem]
    modelFile = open(PREFIX + mfile, 'r')
    (species, logics, kofe, defI, optI) = readModel(modelFile)
    modelFile.close()
    expFile = open(PREFIX + efile, 'r')
    (exps, states) = readExp(expFile)
    expFile.close();

    # convinient variables
    bitlen = len(species)
    kos, fes = kofe['KO'], kofe['FE']
    
    if debug:
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
    for c, s in enumerate(species):
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
    for c, s in enumerate(species):
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
                  L_[s].sort(), isExpOf2(solver.model()[L_[s]]))
            
    # 2. Interactions limit (only for optional interactions)
    allopts = [] # all optional interactions
    if interactions_limit:
        for c, s in enumerate(species):
            actn, repn = map(len,optI[c]) # nums of ats and reps
            # fill the allOpt list with all optional interactions
            if A_[s]:
                a, b = A_[s].size() - actn, A_[s].size()
                allopts.extend([ Extract(i,i,A_[s]) for i in range(a,b) ])
            if R_[s]:
                a, b = R_[s].size() - repn, R_[s].size()
                allopts.extend([ Extract(i,i,R_[s]) for i in range(a, b) ])
        # constraints that all selected nums of inters are less than limit
        solver.add(ULE(sum([ZeroExt(6, bv) for bv in allopts]),
                       interactions_limit))
        if detail:
            print '>> #2 Interactions limit %d ADDED. %s' \
                %(interactions_limit,solver.check())
    else:
        if detail: print '>> #2 Interactions limit: NOT SET.'
        
    # 3. Define Transition ralationship. It is like a macro.
    T = lambda q_old, q_new, ko, fe: \
        And([ And([ Implies(Extract(l, l, L_[s]) == 1, # if logic l is selected
                            (Extract(c,c,q_new)==1) == f_[s][l](q_old, ko, fe))
                    for l in range(L_[s].size())])
              for c, s in enumerate(species) ])
    # this manner is a little bit more quicker
    bunchT = lambda qs, ko, fe: \
        And([ And([ Implies(Extract(l, l, L_[s]) == 1, # if logic l is selected
                            And([ (Extract(c, c, qs[t+1])==1) == \
                                  f_[s][l](qs[t], ko, fe)
                                  for t in range(STEP-1) ]) )
                    for l in range(L_[s].size()) ])
              for c, s in enumerate(species) ])
    if detail: print '>> #3 Transition relation macro: SET.'

    # 4. Experimental constrains
    if detail: print '>> #4 Applying experimental constraints:'
    total_exp = len(exps)
    count_exp = 0
    for exp in exps:
        count_exp += 1
        # build path
        ko_exp = BitVec(exp + '_KO', len(kos) or 1)
        fe_exp = BitVec(exp + '_FE', len(fes) or 1)
        path = [ BitVec(exp + '_%d'%t, bitlen) for t in range(STEP) ]
        #solver.add(*[ T(path[t], path[t+1], KO, FE) for t in range(STEP-1) ])
        solver.add(bunchT(path, ko_exp, fe_exp))
        # add constrains
        for t, conditions in exps[exp]:
            for cond in conditions:
                for s, value in states[cond]:
                    sl = s.split('_')
                    s = sl[-1]
                    if sl[0] == 'KO':
                        c = kos.index(s)
                        solver.add( Extract(c,c,ko_exp) == value)
                    elif sl[0] == 'FE':
                        c = fes.index(s)
                        solver.add( Extract(c,c,fe_exp) == value)
                    else:
                        c = species.index(s)
                        solver.add( Extract(c,c,path[t]) == value )
        if detail:
            print '>> \t %02d/%d %s added...'%(count_exp, total_exp, exp)
    if detail:
        print ">> All Constrains established. (takes %.1f min)" \
            %((time() - startt)/60)

    # Now get the solutions !
    print '>> ' + '-'* 76 # seperator
    solvingt = lastt = time() # just for timing
    print '>> Start solving: %s'%strftime("%d %b %H:%M",localtime(solvingt))
    
    count = 0
    allAR = filter(lambda x: x, list(A_.values()) + list(R_.values()))
    while solver.check() == sat:
        count += 1
        m = solver.model()
        lastt = time() # time for last model
        print ">> Solution %d: (takes %.1f min)"%(count,(time() - lastt)/60)
        if output:
            printModel(m, A_, R_, L_, species, inters, logics,
                       config = True, model = True)
        if count == solutions_limit: break
        # find different solutions (with different selections of interactions)
        solver.add(Or([ b != (m[b] or 0) for b in allAR]))

    endt = time()
    print '>> %d solutions found. (takes %.1f min to end)' \
        %(count, (endt - lastt)/60)
    print '>> End solving: %s.'%strftime("%d %b %H:%M",localtime(endt))
    print '>> ' + '-'* 76
    solving = conv_time(endt - solvingt)
    total = conv_time(endt - startt)
    print '>> Solving duration:\t%s'%solving
    print '>> Total duration:\t%s'%total
    mailMe('Solutions number:\t%d\nTotal duration:\t%s'%(count, total),
           "Computation Fnishied for '%s'."%problem)
    print '>> ' + '-' * 3 + ' Finished. (mail sent.) ' + '-' * 3

if __name__ == '__main__':
    s = Solver()
    main(s,
         problem = 'find_min_inter',
         solutions_limit = 0,
         interactions_limit = 17,
         debug = True, detail = True, output = True)
