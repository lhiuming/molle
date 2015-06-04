from z3 import *
from itertools import combinations
from operator import and_, or_
from time import strftime
from pprint import pprint

def _sorted_inters(inter_list, code):
    ''' Sorts the inter_list = [('from', 'to', 'positive'), ...] into a dict,
     where keyvalue is a couple of number tuples, wich integer codes as keys.

    e.g.
    {1: ( (2, 3), (5, 9) )} : species 1 is activated by 2 and 3, and repressed
    by 5 and 9.
    '''
    d = dict([(c, ([], [])) for c in code.values()]) # initialization
    for i in inter_list:
        f, t = i[:2]
        idx = (0, 1)[ i[2]=='negative' ]
        d.setdefault(code[t], ([], []))[idx].append(code[f])
    return d

def readModel(f):
    ''' Take a file Object as input, return a tuple of 6 objects:

    species: a list of gene name.
    code   : a dict of gene_name : integer_code. Consistant with species list.
    logics : a dict. { gene_name: list_of_allowed_logic_numbers }
    kofe   : a dict. { "FE": list_of_FEable_gene, "KO": list_of_KOable_gene }
    defI   : a dict of defined interations. Processed by _sorted_inters()
    optI   : a dict of optional interactions.
    '''
    species = []
    code = {}
    logics = {}
    kofe = {'KO':[], 'FE':[]}
    def_inters_list = []
    opt_inters_list = []
  
    # read the components line
    for c in f.readline().strip().split(','):
        if '(' in c :
            gene_ = c[:c.index('(')].strip() # name is before '('
        else:
            gene_ = c.strip() # no logics specified
        gene = filter(lambda x: not x in '+-', gene_)
        mark = filter(lambda x: x in '+-', gene_)
        # add to kofe if the gene has mark
        if('+' in mark): kofe['FE'].append(gene)
        if('-' in mark): kofe['KO'].append(gene)
        # record allows logics
        if '(' in c:
            rules = tuple([int(i)\
                           for i in c[c.index('(')+1:c.index(')')].split()])
        else: # no specific logics
            if gene in ('MEKERK', 'Tcf3'): rules = (16, 17)
            else: rules = tuple(range(16))
        logics[gene] = rules # record the list of allowed functions
        species.append(gene)
    code = dict([ (s, n) for (n, s) in enumerate(species) ])

    # read the interaction lines
    for l in f.readlines(): # loop every line
        l = l.strip().split()
        if(not l): continue # skip empty line
        if(l[-1] == 'optional'): opt_inters_list.append(tuple(l[:3]))
        else: def_inters_list.append(tuple(l[:3]))
    defI = _sorted_inters(def_inters_list, code)
    optI = _sorted_inters(opt_inters_list, code)

    return (species, code, logics, kofe, defI, optI)

# kept from old version
def _addExp(d, name, time_point, state_names_list):
    d.setdefault(name, []).append( (int(time_point), state_names_list) )
  
# kept from old version
def _addState(d, state_name, gene, value):
    d.setdefault(state_name, []).append( (gene, int(value)) )

# kept from old version
def readExp(f):
  '''
  Take the file for experiment constrains, return two dicts:

  exps:   the Experimental constrains for every experiment
  states: records the mapping of shortcut name to node states
  '''
  exps = dict()
  states = dict()

  shortcut = ''
  for l in f.readlines():
    l = l.strip();
    if(not l): continue; # skip empty line

    try: l = l[:l.index('"')] # remove commment
    except ValueError: None
    try: l = l[:l.index(';')] # remove ;
    except ValueError: None

    if(shortcut): # inside the braket { }
      if(l[0] == '{'): continue # skip left bracket
      elif(l[0] == '}'): shortcut = '' # exit the braket;
      else:
        (left, right) = l.split('=');
        name = left.strip();
        value = right.split()[0];
        _addState(states, shortcut, name, value); # record configuration
    l = l.split();
    if(l[0] == "//"): continue # comment line
    elif(l[0] == "under"): _addExp(exps, l[1], l[3], l[4:]) # recordexp
    elif(l[0] == "let"):
     shortcut = l[1]; # ready to enter the braket
     try: shortcut = shortcut[:shortcut.index(':')]
     except ValueError: None

  return (exps, states);

def _sublist_generator(l, high=None, low=0):
    ''' Yields a sublist. '''
    if not high: high = len(l)
    for r in range(low, high + 1):
        g = combinations(l, r)
        for sl in g: yield sl

def generateInterCombi(defs, opts):
    ''' Yielf one of *all* iteraction combinations, represented by a couple
    of tuples, containing activators and represors.
    '''
    defAct, defRep = defs
    optAct, optRep = opts
    for act in _sublist_generator(optAct):
        for rep in _sublist_generator(optRep):
            yield ( tuple(defAct + list(act)), tuple(defRep + list(rep)) )

# kept from older version
def _compati_func(acrp, funclist):
  '''
  Return a filtered tuple of function numers. Only meaningful and allowed
  function for the node is included, accroding the number of activators and
  repressors. acrp is the tuples of activators list and repressors list.
  funclist is a tuple for allowed functions.
  '''
  acn = len(acrp[0]) # number of activators
  rpn = len(acrp[1]) # number of repressors

  if(acn == 0): # when no activator
    if(rpn == 0): return (-1, ) # -1 means always False
    if(rpn == 1): return filter(lambda x: x == 17, funclist) or (-1, )
    if(rpn >= 2): return filter(lambda x: x >= 16, funclist) or (-1, )
  if(acn == 1): # when only one activator. assumes funclist = 0 ~ 15
    if(rpn == 0): return (1, ) # 1 means only activated when activator presents
    if(rpn == 1): return (0, 2, 8) # this is representative
    if(rpn >= 2): return (0, 2, 4, 8) # representative
  if(acn >= 2): # when two activators. assumes funclist = 0 ~ 15
    if(rpn == 0): return (0, 1) # representative
    if(rpn == 1): return (0, 2, 4, 8, 10, 14) # representative
    if(rpn >= 2): return (0, 2, 4, 6, 8, 10, 12, 14) # representative

zero = BitVecVal(0, 1)
one = BitVecVal(1, 1)

def any(bvs):
    return reduce(or_, bvs, 0)

def all(bvs):
    return reduce(and_, bvs, 1)
    
def _concat(bvs):
    if len(bvs) == 1: return bvs[0]
    else: return Concat(bvs)
    
def _create_bit_rule(num, act, rep, A, R):
    ''' Create the update rule that return bit-vector of length 1. '''
    if act:
        act = _concat(act)
        aa = act & A
    if rep:
        rep = _concat(rep)
        rr = rep & R
    
    if act:
        if not rep: # not repressor, but have activators
            if num%2: return act == A
            else: return aa != 0
        else: # both activators and repressors present
            if num < 2: return False
            elif num < 4: return And(aa != 0, rr == 0)
            elif num<6: return And(act == A, rep != R)
            elif num<8: return And(rr != 0, rep != R)
            elif num<10: return act == A
            elif num<12: return Or(act == A, And(aa != 0, rr == 0))
            elif num<14: return Or(act == A, And(aa != 0, rep != R))
            elif num<16: return aa != 0
            else: return False
    if rep: # no activator but have repressors
        if num==16: return And(rr != 0, rep != R)
        elif num==17: return rr == 0
        else: return False
    return False

# kept from older version
def _with_KOFE(node, kofe, prec):
  if(node in kofe['KO']):
    if(node in kofe['FE']):
      return lambda x: Or(prec['FE'][node], And(Not(prec['KO'][node]), x))
    else: return lambda x: And(Not(prec['KO'][node]), x) # only ko
  if(node in kofe['FE']):
    return lambda x: Or(prec['FE'][node], x)
  else: return lambda x: x; # no kofe

def makeFunction(acts, reps, logic, A, R):
    ''' Makes a function that takes q, A, R, and return a coresponding z3 expr.
    A is the acticators-selecting bit-vector, R for repressors.
    '''
    return lambda q: _create_bit_rule(logic,
                                      [Extract(i, i, q) for i in acts],
                                      [Extract(i, i, q) for i in reps],
                                      A, R)


### Output Utilities ###
#########################
boolf = BoolVal(False)

def _Or(l):
    if(not l): return boolf
    if(len(l) == 1): return l[0]
    else: return Or(l);

def _And(l):
    if(not l): return boolf
    if(len(l) == 1): return l[0]
    else: return And(l);

def _create_sym_rule(num, act, rep):
    if(num == -1): return boolf
    if(num < 2 and rep): return boolf
    if(num > 15 and act): return boolf
    
    actt = [Bool(node) for node in act]
    rept = [Bool(node) for node in rep]
        
    if(num > 1 and not rep): return (_And, _Or)[num % 2](actt)

    if(num == 0): return _And(actt)
    elif(num == 1): return _Or(actt)
    elif(num < 4): return And(_Or(actt), Not(_Or(rept)))
    elif(num < 6): return And(_And(actt), Not(_And(rept)));
    elif(num < 8): return And(_Or(actt), Not(_And(rept)))
    elif(num < 10): return _And(actt)
    elif(num < 12): return Or(_And(actt), And(_Or(actt), Not(_Or(rept))))
    elif(num < 14): return Or(_And(actt), And(_Or(actt), Not(_And(rept))))
    elif(num < 16): return _Or(actt)
    elif(num == 16): return And(_Or(rept), Not(_And(rept)))
    elif(num == 17): return Not(_Or(rept));

def checkBit(i, bv):
    return simplify(Extract(i, i, bv)).as_long()

def printModel(m, A_, R_, L_, species, code, inters, config = True):
    ''' Print the solved model nicely. '''
    # getting model details
    A = {}; R = {}; L = {}
    for s in species:
        c = code[s]
        L[s] = len(bin(m[L_[s]].as_long()).lstrip('0b')) - 1
        if A_[s]:
            actB = m[A_[s]]; l = actB.size() - 1
            A[s] = [species[n] for i, n in enumerate(inters[c][0]) \
                    if checkBit(l-i, actB) ]
        else: A[s] = []
        if R_[s]:
            repB = m[R_[s]]
            if not repB: print m
            l = repB.size() - 1
            R[s] = [species[n] for i, n in enumerate(inters[c][1]) \
                    if checkBit(l-i, repB) ]
        else: R[s] = []
    # printing details
    print '>>'
    if config:
        print ">>\tConfigurations: "
        for s in species:
            print ">>\t\t%s:%d%s%s" \
            %(s, L[s],
              A[s] and '\t<- ' + ','.join(A[s]) or '',
              R[s] and '\t|- ' + ','.join(R[s]) or '')
    print ">>\tModel: "
    for s in species: print ">>\t\t%s' = %s" \
        %(s,simplify( _create_sym_rule(L[s], A[s], R[s])))


### Debugging Secntions ###
###########################
if __name__ == "__main__":
  if __debug__:
      print ">> testing file reading: "
      modelFile = open('examplefiles/SimpleFourComponentModel.txt', 'r')
      expFile = open('examplefiles/CertainInteractionRequired.txt', 'r')
      (species, code, logics, kofe, defI, optI) = readModel(modelFile)
      (exps, states) = readExp(expFile)
      modelFile.close(); expFile.close()
      print "speceis code: ", code
      print "optional interactions: ", optI

      print ">> testing getInterCombi():"
      testDef = (['a1', 'a2'], ['r1', 'r2'])
      testOpt = (['aaa'], [])
      pprint([i for i in generateInterCombi(testDef, testOpt)])

      print ">> testing makeFunction():"
      testq = BitVec('testq', 8)
      testI = ( (1, 3), (6, 7) )
      print makeFunction(testI, 1)(testq)
      print makeFunction(testI, 16)(testq)
      print simplify(makeFunction(testI, 7)(testq) == 1)


