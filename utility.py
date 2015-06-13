from z3 import *
from operator import or_
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

zero = BitVecVal(0, 1)

def Any(bvs):
    return reduce(or_, bvs, zero)

def _concat(bvs):
    if len(bvs) == 1: return bvs[0]
    else: return Concat(bvs)
    
def _create_bit_rule(num, act, rep, A, R):
    ''' Create the update rule that return bit-vector of length 1. '''
    # initialization
    if act:
        act = _concat(act)
        aa = act & A
    if rep:
        rep = _concat(rep)
        rr = rep & R
    # creating result
    if act:
        if not rep: # not repressor, but have activators
            if num%2 == 0: return act == A
            else: return aa != 0
        else: # both activators and repressors present
            if num == 0: return And(R == 0, act == A)
            elif num==1: return And(R == 0, aa != 0)
            elif num==2: return Or( And(R == 0, act == A),
                                    And(R != 0, aa != 0, rr == 0) )
            elif num==3: return Or( And(R == 0, aa != 0),
                                    And(R != 0,aa != 0, rr == 0) )
            elif num==4: return Or( And(R == 0, act == A),
                                    And(R != 0, act == A, rep != R) )
            elif num==5: return Or( And(R == 0, aa != 0),
                                    And(R != 0, act == A, rep != R) )
            elif num==6: return Or( And(R == 0, act == A),
                                    And(R != 0, aa != 0, rep != R) )
            elif num==7: return Or( And(R == 0, aa != 0),
                                    And(R != 0, aa != 0, rep != R) )
            elif num==8: return Or( And(R == 0, act == A),
                                    And(R != 0, act == A) )
            elif num==9: return Or( And(R == 0, aa != 0),
                                    And(R != 0, act == A) )
            elif num==10: return Or( And(R == 0, act == A),
                                     And(R != 0,
                                         Or(act == A, And(aa != 0, rr == 0))))
            elif num==11: return Or( And(R == 0, aa != 0),
                                     And(R != 0,
                                         Or(act == A, And(aa != 0, rr == 0))))
            elif num==12: return Or( And(R == 0, act == A),
                                     And(R != 0,
                                         Or(act == A, And(aa != 0, rep != R))))
            elif num==13: return Or( And(R == 0, aa != 0),
                                     And(R != 0,
                                         Or(act == A, And(aa != 0, rep != R))))
            elif num==14: return Or( And(R == 0, act == A),
                                     And(R != 0, aa != 0) )
            elif num==15: return aa != 0
            else: return False
    if rep: # no activator but have repressors
        if num==16: return And(rr != 0, rep != R)
        elif num==17: return rr == 0
        else: return False
    return False

def _with_kofe(kofe_idx, ko, fe, expr):
    koc, fec = kofe_idx
    if koc:
        ko = Extract(koc-1,koc-1,ko) == 1 # a trick to avoid 0 == False
        if fec:
            fe = Extract(fec-1,fec-1,fe) == 1
            return Or(fe, And(Not(ko), expr))
        else: return And(Not(ko), expr)
    elif fec:
        fe = Extract(fec-1,fec-1,fe) == 1
        return Or(fe, expr)
    else: return expr

def makeFunction(acts, reps, kofe_index, logic, A, R):
    ''' Makes a function that takes q, A, R, and return a coresponding z3 expr.
    A is the acticators-selecting bit-vector, R for repressors.
    '''
    return lambda q, ko, fe: \
        _with_kofe(kofe_index, ko, fe,
                   _create_bit_rule(logic, [Extract(i,i,q) for i in acts],
                                    [Extract(i, i, q) for i in reps], A, R))

def isExpOf2(bvv):
    return len(filter(lambda x: x == '1', bin(bvv.as_long()))) == 1

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
    if num < 0: return Bool('Strange')
    if act:
        actt = [Bool(node) for node in act]
    if rep:
        rept = [Bool(node) for node in rep]

    if act:
        if not rep:
            if num%2 == 0: return _And(actt)
            else: return _Or(actt)
        elif num == 0: return boolf
        elif num == 1: return boolf
        elif(num < 4): return And(_Or(actt), Not(_Or(rept)))
        elif(num < 6): return And(_And(actt), Not(_And(rept)));
        elif(num < 8): return And(_Or(actt), Not(_And(rept)))
        elif(num < 10): return _And(actt)
        elif(num < 12): return Or(_And(actt), And(_Or(actt), Not(_Or(rept))))
        elif(num < 14): return Or(_And(actt), And(_Or(actt), Not(_And(rept))))
        elif(num < 16): return _Or(actt)
        else: return boolf
    if rep:
        if num == 16: return And(_Or(rept), Not(_And(rept)))
        elif num==17: return Not(_Or(rept));
        else: return boolf
    return boolf

def checkBit(i, bv):
    return simplify(Extract(i, i, bv)).as_long()

def bv2logic(lbvv, llist):
    ''' convert a bit-vector to a integer, as logic function number.'''
    assert isExpOf2(lbvv)
    lcode = len(bin(lbvv.as_long()).lstrip('0b')) - 1
    return llist[lcode]

def printModel(m, A_, R_, L_, species, code, inters, logics,
               config = True, model = True):
    ''' Print the solved model nicely. '''
    # getting model details
    A = {}; R = {}; L = {}
    for s in species:
        c = code[s]
        L[s] = bv2logic(m[L_[s]], logics[s])
        if A_[s]:
            actB = m[A_[s]]; l = actB.size() - 1
            A[s] = [species[n] for i, n in enumerate(inters[c][0]) \
                    if checkBit(l-i, actB) ]
        else: A[s] = []
        if R_[s]:
            repB = m[R_[s]]; l = repB.size() - 1
            R[s] = [species[n] for i, n in enumerate(inters[c][1]) \
                    if checkBit(l-i, repB) ]
        else: R[s] = []
    # printing the model
    print '>>'
    if config:
        print ">>\tConfigurations: "
        for s in species:
            print ">>\t\t%s:%d%s%s" \
            %(s, L[s],
              A[s] and '\t<- ' + ','.join(A[s]) or '',
              R[s] and '\t|- ' + ','.join(R[s]) or '')
    if model:
        print ">>\tModel: "
        for s in species: print ">>\t\t%s' = %s" \
            %(s,simplify( _create_sym_rule(L[s], A[s], R[s]) ))
    print '>>'


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


