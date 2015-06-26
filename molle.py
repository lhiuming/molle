from z3 import *
from utility import *
from pprint import pprint


STEP = 20 # trajectory length of the model

class ABN:
    ''' The object receive input files: model = Model(mfile, efile) .

    build(): Build the model by model.build(interaction_limit).
    solve(): Get the solutions by model.solve(solutions_limit) after build.
    '''
    
    def __init__(self, mfile, efile):
        (self.species, self.logics, self.kofe,
         self.defI, self.optI) = readModel(mfile)
        (self.exps, self.states) = readExp(efile)

    def build(self, ilimit=0, detail=True, debug=False):
        ''' Add all constrains, and set self.solver. '''
        bitlen = len(self.species)
        kos, fes = self.kofe['KO'], self.kofe['FE']        
        solver = Solver()
    
        if debug:
            print '>> Defined Iteractions:'; pprint(self.defI)
            print '>> Optional Interactions:'; pprint(self.optI)
            print '>> KO: ', kos
            print '>> FE: ', fes

        # 0. Encoding functions
        self.A_ = {} # of activator/activating-interaction selection BitVec
        self.R_ = {} # of repressor/repressing-interaction selection BitVec
        self.L_ = {} # of Logic-selecting BitVec for speciew
        self.f_ = {} # devices/modules/functions
        self.inters = {} # record all avalible acts and reps for every specie
        for c, s in enumerate(self.species):
            acts = self.optI[c][0] + self.defI[c][0] # Concat is from L to R
            reps = self.optI[c][1] + self.defI[c][1]
            self.inters[c] = (acts, reps)
            # filter the useless functiono
            self.logics[s] = compati(self.logics[s], len(acts), len(reps))
            
            # creating Act and Rep selecting BitVec
            if acts: self.A_[s] = BitVec('Act_' + s, len(acts))
            else: self.A_[s] = None
            if reps: self.R_[s] = BitVec('Rep_' + s, len(reps))
            else: self.R_[s] = None
            
            # create logic-selecting BitVec
            self.L_[s] = BitVec('Logic_' + s, len(self.logics[s]))
            
            # make the functions
            kofe_index = (s in kos and kos.index(s) + 1,
                          s in fes and fes.index(s) + 1)
            self.f_[s] = [ makeFunction(acts, reps, kofe_index,
                                        l, self.A_[s], self.R_[s])
                           for l in self.logics[s] ]

        # 1. Modeling Constrains
        if detail:
            print '>> #1 Adding modeling constrains: '
        for c, s in enumerate(self.species):
            # INTER: defined activators and repressors must be selected
            defactn, defrepn = map(len, self.defI[c])
            solver.add([ 1== Extract(i,i,self.A_[s]) for i in range(defactn) ])
            solver.add([ 1== Extract(i,i,self.R_[s]) for i in range(defrepn) ])
            
            # LOGIC: only one logic is selected
            logic_i = range(self.L_[s].size())
            solver.add(
                0 == Any([ Extract(i,i,self.L_[s]) & Extract(j,j,self.L_[s])
                           for i in logic_i for j in logic_i if i != j]) )
            # LOGIC: must select one logic
            solver.add(self.L_[s] != 0)
            
            if debug:
                assert solver.check() == sat
                tmp = lambda b: b and b.sort() or None
                print '>> \tConstraints %s:\tAct=%s,\tRep=%s,\tLog=%s(%s).' \
                    %(s, tmp(self.A_[s]), tmp(self.R_[s]),
                      self.L_[s].sort(), isExpOf2(solver.model()[self.L_[s]]))


        # 2. Interactions limit (only for optional interactions)
        opts = [] # all optional interactions
        if ilimit:
            for c, s in enumerate(self.species):
                actn, repn = map(len,self.optI[c]) # nums of ats and reps
                # fill the allOpt list with all optional interactions
                if self.A_[s]:
                    a, b = self.A_[s].size() - actn, self.A_[s].size()
                    opts.extend([ Extract(i,i,self.A_[s]) for i in range(a,b)])
                if self.R_[s]:
                    a, b = self.R_[s].size() - repn, self.R_[s].size()
                    opts.extend([ Extract(i,i,self.R_[s]) for i in range(a,b)])
                # all selected nums of inters are less than limit
                solver.add(ULE(sum([ZeroExt(6, bv) for bv in allopts]),
                               ilimit))
            if detail:
                print '>> #2 Interactions limit %d ADDED. %s' \
                    %(ilimit, solver.check())
        else:
            if detail: print '>> #2 Interactions limit: NOT SET.'
        
        # 3. Experimental constrains
        if detail:
            print '>> #3 Applying experimental constraints:'
        total_exp = len(self.exps)
        count_exp = 0
        for name, exp in self.exps.items():
            count_exp += 1
            # build updating path
            ko_exp = BitVec(name + '_KO', len(kos) or 1)
            fe_exp = BitVec(name + '_FE', len(fes) or 1)
            path = [ BitVec(name + '_%d'%t, bitlen) for t in range(STEP) ]
            #solver.add(*[ T(path[t], path[t+1], KO, FE)
            #              for t in range(STEP-1) ])
            solver.add(self.bunchT(path, ko_exp, fe_exp))
            # add constrains
            for t, conditions in exp:
                for cond in conditions:
                    for s, value in self.states[cond]:
                        if s[:3] == 'KO_':
                            c = kos.index(s[3:])
                            solver.add( Extract(c,c,ko_exp) == value )
                        elif s[:3] == 'FE_':
                            c = fes.index(s[3:])
                            solver.add( Extract(c,c,fe_exp) == value )
                        else:
                            c = self.species.index(s)
                            solver.add( Extract(c,c,path[t]) == value )
            if detail:
                print '>> \t %02d/%d %s added...'%(count_exp, total_exp, name)

        self.solver = solver

    def T(self, q_old, q_new, ko, fe):
        ''' Define Transition ralationship. It is like a macro.'''
        return \
          And([
            And([Implies(Extract(l, l, self.L_[s]) == 1,
                         (Extract(c,c,q_new)==1) == self.f_[s][l](q_old,ko,fe))
                 for l in range(self.L_[s].size())])
            for c, s in enumerate(self.species) ])

    def bunchT(self, qs, ko, fe):
        ''' this manner is a little bit more quicker.'''
        return \
          And([
            And([ Implies(Extract(l, l, self.L_[s]) == 1,
                          And([ (Extract(c, c, qs[t+1])==1) == \
                                self.f_[s][l](qs[t], ko, fe)
                                for t in range(STEP-1)]) )
                  for l in range(self.L_[s].size()) ])
            for c, s in enumerate(self.species) ])

    def solve(self):
        ''' Get the solutions. return an iterator.'''
        allAR = filter(None, list(self.A_.values()) + list(self.R_.values()))
        while self.solver.check() == sat:
            m = self.solver.model()
            yield Solution(m, self.A_, self.R_, self.L_,
                           self.species, self.inters, self.logics)
            self.solver.add(Or([ b != (m[b] or 0) for b in allAR]))

class Solution:
    ''' A printable solution object. '''
    def __init__(self, m, A_, R_, L_, species, itrs, lgcs):
        self.m = m
        self.species = species
        self.A, self.R, self.L = getDetail(m, A_, R_, L_, species, itrs, lgcs)

    def output(self, config=True, model=True):
        printModel(self.species, self.A, self.R, self.L, config, model)
        
        
    
