from z3 import *
from itertools import combinations, product
from time import strftime
from pprint import pprint

def readModel(f):
  '''
  Take a file Object as input, return a 4 objects:
    comps    : a dict. { gene_name: list_of_allowed_logic_function_number }
    kofe     : a dict. { "FE": list_of_FEable_gene, "KO": list_of_KOable_gene }
    defInters: a list of defined interations. [ ("from", "to", "positive" ), ]
    optInters: a list of optional interactions.
  '''
  comps = dict()
  kofe = {'KO':[], 'FE':[]}
  defInters = []
  optInters = []
  
  # read the components line
  for c in f.readline().strip().split(','):
    if('(' in c): # functions are specified
      gene_ = c[:c.index('(')].strip() # name is befor '('
    else: # no '('  -> no specific functions
      gene_ = c.strip()
      
    gene = filter(lambda x: not x in '+-', gene_) # get the gene name
    mark = filter(lambda x: x in '+-', gene_) # get the +- mark
    for m in mark: # add to kofe if the gene has mark
      if(m == '+'): kofe['FE'].append(gene)
      elif(m == '-'): kofe['KO'].append(gene)

    if(not '(' in c): # no specific funcions
      if(gene in ('MEKERK', 'Tcf3')): rules = (16, 17)
      else: rules = tuple(range(16))
    else: # functions specified, then read them in the '(' ')'
      rules = tuple([int(i) for i in c[c.index('(')+1:c.index(')')].split()])

    comps[gene] = rules # record the list of allowed functions

  # read the interaction lines
  for l in f.readlines(): # loop every line
    l = l.strip().split()
    if(not l): continue # skip empty line
    if(l[-1] == 'optional'): optInters.append(tuple(l[:3]))
    else: defInters.append(tuple(l[:3]))

  return (comps, kofe, defInters, optInters)

def _addExp(d, name, time_point, state_names_list):
  _appendEntry(d, name, (int(time_point), state_names_list))

def _addState(d, state_name, gene, value):
  _appendEntry(d, state_name, (gene, int(value)))

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

def preCon(comps, kofe = None, step = 20):
  d = dict();
  for node in comps:
    d[node] = [Bool(node + '_' + str(t)) for t in range(step)];
  if(kofe):
    d['KO'] = dict([(node, Bool('KO_' + node)) for node in kofe["KO"]])
    d['FE'] = dict([(node, Bool('FE_' + node)) for node in kofe["FE"]])
  return d;

def _get_sublist(l, limit):
  '''
  Generator for non-empty sublists of list l. It will generate sublists in
  a decreasing size, first tuple(l), last a empty tuple.
  '''
  if(limit): rg = range( min(len(l), limit) + 1 )
  else: rg = range(len(l) + 1)
  rg.reverse() # return the longest sublist first
  for r in rg:
    for combi in combinations(l, r):
      yield combi;

def _appendEntry(d, key, entry):
  try:
    d[key].append(entry)
  except KeyError:
    d[key] = [entry]

def _mergeInters(d, bemerged):
  for node in bemerged:
    try: d[node].extend(bemerged[node]);
    except KeyError:
      d[node] = bemerged[node][:]; # a shallow copy;

def _intoInters(interList):
  d = dict()
  for l in interList:
    _appendEntry(d, l[1], ( l[0], (l[2] == "positive") and 1 or -1 ))
  return d

def _construct_graph(l, defI = None, allNodes = None):
  '''
  Take a list of directed and signed interactions, return a dictionary of nodes
  with inputs list as value. If defInters are provide, it will be merge in the
  resulting graph.
  '''
  d = _intoInters(l)
  if(defI):
    _mergeInters(d, defI); # records in defI will be merged in d
  if(allNodes):
    for node in allNodes:
      if(node not in d): d[node] = [];
  return d;

def getGraph(comps, optI, defI, interLimit):
  '''
  Generator. Yield a subpgraph from all interactions. Defined
  interactions are bound to be present.
  '''
  complist = comps.keys();
  defInters = _intoInters(defI)

  for sl in _get_sublist(optI, max(interLimit, 0)): # get a sub-list of inters
      yield _construct_graph(sl, defInters, complist); # return a sub-graph

def sortGraph(graph):
  '''
  Reconstruct the graph dict in to a more compact dict.
  '''
  d = dict();
  for node in graph:
    act = filter(lambda x: x[1] > 0, graph[node]);
    rep = filter(lambda x: x[1] < 0, graph[node]);
    d[node] = ([x[0] for x in act], [x[0] for x in rep]);
  return d;

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

def getFunction(comps, sgraph, compact = False):
  '''
  Generator a combination of all possible combinations of allowed function for
  every node. This generator will also check the compatability of
  rules and input nodes. The yielding is a dictionary, with node:function.
  '''
  nodes = comps.keys()
  rules = dict()
  i = 1
  for n in nodes:
    if(compact) :rules[n] = list(_compati_func(sgraph[n], comps[n]))
    else: rules[n] = list(comps[n])
    rules[n].reverse()
    i *= len(rules[n]);
  yield (i, rules) # return the configuration as the first yield
  for c in product(*[rules[node] for node in nodes]):
    yield dict(zip(nodes, c));

def _Or(l):
  if(not l): return False
  if(len(l) == 1): return l[0]
  else: return Or(l);

def _And(l):
  if(not l): return False
  if(len(l) == 1): return l[0]
  else: return And(l);

def _create_rule(num, act, rep, precon, t = None):
  if(num == -1): return False
  if(num < 2 and rep): return False
  if(num > 15 and act): return False

  actt = [(precon[node])[t] for node in act]
  rept = [(precon[node])[t] for node in rep]
  
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

def _with_KOFE(node, kofe, prec):
  if(node in kofe['KO']):
    if(node in kofe['FE']):
      return lambda x: Or(prec['FE'][node], And(Not(prec['KO'][node]), x))
    else: return lambda x: And(Not(prec['KO'][node]), x) # only ko
  if(node in kofe['FE']):
    return lambda x: Or(prec['FE'][node], x)
  else: return lambda x: x; # no kofe

def _add_function(s, fnum, node, kofe, acrp, prec = None, step = 20):
  (ac, rp) = acrp;
  tune = _with_KOFE(node, kofe, prec);
  for t in range(1, step): # 19 transitions from t-1 to t
    s.add(Bool(node + '_' + str(t)) == tune(
               _create_rule(fnum, ac, rp, prec, t = t-1)));

def applyFunctions(s, funcd, kofe, sgraph, prec = None, step = 20):
  '''
  Take the solver, function number, and a sorted graph.
  '''
  for node in funcd:
    _add_function(s, funcd[node], node, kofe, sgraph[node], prec, step);
    
def addConstrains(s, exp, states, precon):
  '''
  exp is a list of tuple, in the form (t, [statenames]). states is the dict
  for shortcut of conditions.
  '''
  for (t, shortcuts) in exp:
    for sc in shortcuts:
      for (node, value) in states[sc]:
        par = node.split('_')
        if(par[0] in ('KO', 'FE')):
          s.add(precon[par[0]][par[1]] == bool(value))
        else:
          s.add(precon[node][t] == bool(value));
        

### Output Utilities ###
#########################

def _reprModel(sgraph, funcd):
  '''
  Return a readable presentation of model for output. Resemble those ouput by
  Microsoft's webapp.
  '''
  s = Solver()
  for node in sgraph:
    (ac, rp) = sgraph[node]
    s.add(Bool(node + "'") == _create_rule(funcd[node], ac, rp))
  return s

def printModel(sgraph, funcd):
  ''' Useful for print model on stdout. '''
  pprint(_reprModel(sgraph, funcd))
  print funcd, '\n'

def outputModel(sgraph, funcd, fileName, count):
  '''
  Record the model in a file specified by fileName.
  '''
  s = _reprModel(sgraph, funcd)
  with open(fileName, "a") as f:
    f.write(strftime(">> %d %b %H:%M:%S  "))
    if(count): f.write("The %dth solution is: \n" %count)
    f.write(str(s) + '\n')
    f.write(">> Interactoins and functions are: \n")
    f.write(str(sgraph) + '\n')
    f.write(str(funcd) + '\n\n')
    f.close()


