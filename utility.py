from z3 import *
from itertools import combinations
from operator import and_
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
    for i in inters:
        f, t = i[:2]
        idx = (0, 1)[ i[2]=='negative' ]
        d.setdefault(code[t], ([], []))[idx].append(code[f])
    return d

def readModel(f):
    ''' Take a file Object as input, return 6 objects:

    species: a list of gene name.
    code   : a dict of gene_name : integer_code. consistant with species list.
    logics : a dict. { gene_name: list_of_allowed_logic_function_number }
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

def _addExp(d, name, time_point, state_names_list):
  _appendEntry(d, name, (int(time_point), state_names_list))

def _addState(d, state_name, gene, value):
  _appendEntry(d, state_name, (gene, int(value)))

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

def getInterCombi(defs, opts):
    ''' Return a list of *all* iteraction combinations, represented by a couple
    of tuples, containing activators and represors.
    '''
    defAct, defRep = defs
    optAct, optRep = opts
    return [(tuple(defAct + list(act)), tuple(defRep + list(rep))) \
            for act in _sublist_generator(optAct) \
            for rep in _sublist_generator(optRep)]

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

# kept from older version
def _Or(l):
  if(not l): return False
  if(len(l) == 1): return l[0]
  else: return Or(l);

def _And(l):
  if(not l): return False
  if(len(l) == 1): return l[0]
  else: return And(l);

def _concat(l):
    if not l:return BitVecVal(0, 1)
    if len(l) == 1: return l[0]
    return Concat(l)
  
def _create_rule(num, act, rep):
    ''' Create the update rule for making functions. '''
    if num == -1: return False # -1 is the reserved num
    if num < 2 and rep: return False
    if num > 15 and act: return False
    
    actB = _concat(act)
    repB = _concat(rep)

    if num==0: return ~actB == 0
    elif num==1: return actB != 0
    elif num<4: return And(actB != 0, repB == 0)
    elif num<6: return And(~actB == 0, ~repB != 0)
    elif num<8: return And(actB != 0, ~repB != 0)
    elif num<10: return ~actB == 0
    elif num<12: return Or(~actB == 0, And(actB != 0, repB == 0))
    elif num<14: return Or(~actB == 0, And(actB != 0, ~repB != 0))
    elif num<16: return actB != 0
    elif num==16: return And(repB != 0, ~repB != 0)
    elif num==17: return repB == 0
    
# kept from older version
def _with_KOFE(node, kofe, prec):
  if(node in kofe['KO']):
    if(node in kofe['FE']):
      return lambda x: Or(prec['FE'][node], And(Not(prec['KO'][node]), x))
    else: return lambda x: And(Not(prec['KO'][node]), x) # only ko
  if(node in kofe['FE']):
    return lambda x: Or(prec['FE'][node], x)
  else: return lambda x: x; # no kofe

def makeFunction(inter, logic_num):
    ''' make a function that takes a q, and return a coresponding z3 expr.'''
    return lambda q: _create_rule(logic_num,
                                  [Extract(i, i, q) for i in inter[0]],
                                  [Extract(i, i, q) for i in inter[1]])
  
### Output Utilities ###
#########################

def printModel():
    return None


### Debugging Secntions ###
###########################
if __name__ == "__main__":
  if __debug__:
    print ">> testing sortedInters():"
    print sortedInters([['a', 'b', 'positive'], ['b', 'c', 'negative']],
                       {'a':0, 'b':1, 'c':2})

    print ">> testing getInterCombi():"
    testDef = (['a1', 'a2'], ['r1', 'r2'])
    testOpt = (['aaa'], [])
    print getInterCombi(testDef, testOpt)

    print ">> testing makeFunction():"
    testq = BitVec('testq', 8)
    testI = ( (1, 3), (6, 7) )
    print makeFunction(testI, 1)(testq)
    print makeFunction(testI, 16)(testq)
    print makeFunction(testI, 7)(testq)


