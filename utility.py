from z3 import *
from itertools import combinations
from itertools import product

def __init__():
  return;

def appendEntry(d, k, entry):
  try: d[k].append(entry);
  except KeyError:
    d[k] = [];
    appendEntry(d, k, entry);

def mergeInters(d, bemerged):
  for node in bemerged:
    try: d[node].extend(bemerged[node]);
    except KeyError:
      d[node] = bemerged[node][:]; # a shallow copy;

def addInter(d, a, b, p):
  appendEntry(d, b, (a,(p == 'positive') and 1 or -1));

def readModel(f):
  '''
  Take a file as input, return three dicts, which are allowed rules for 
  components,  defined interactions and optional interactions, seperately.
  The entry of interactions dict collect "from" nodes.
  '''
  comps = dict();
  kofe = {'KO':[], 'FE':[]};
  defInters = dict();
  optInters = dict();

  for c in f.readline().split(','): # read the components line
    c = c.strip();
    gene = c[:c.index('(')]; # the first word is the name
    if(gene[-1] == '+'): 
      gene = gene[:-1];
      kofe['FE'].append(gene.strip('-'));
    if(gene[-1] == '-'):
      gene = gene[:-1];
      kofe['KO'].append(gene.strip('+'));
    rulelist = c[c.index('(')+1:c.index(')')];
    comps[gene] = [int(i) for i in rulelist.split()]; # store allowed rules num
    
  for l in f.readlines(): # loop every line for interaction
    l = l.strip().split();
    if(not l): continue;
    if(l[-1] == 'optional'): addInter(optInters, l[0], l[1], l[2])
    else: addInter(defInters, l[0], l[1], l[2]); # defined interaction

  return (comps, kofe, defInters, optInters);

def addExp(d, name, t, stateList):
  appendEntry(d, name, (t, stateList));

def addState(d, name, c, value):
  appendEntry(d, name, (c, value));

def readExp(f):
  '''
  Take the file for experiment constrains. Return two dicts, exps and states.
  exps is the Experimental constrains for every experiment; states records
  the mapping of shortcut name to node states.
  '''
  exps = dict();
  states = dict();

  shortcut = '';
  for l in f.readlines():
    l = l.strip();
    if(not l): continue; # skip empty line

    try: l = l[:l.index('"')] # remove commment
    except ValueError: None;
    try: l = l[:l.index(';')] # remove ;
    except ValueError: None;

    if(shortcut): # inside the braket { }
      if(l[0] == '{'): continue # skip left bracket
      elif(l[0] == '}'): shortcut = '' # exit the braket;
      else:
        (left, right) = l.split('=');
        name = left.strip();
        value = right.split()[0];
        addState(states, shortcut, name, value); # record configuration
    l = l.split();
    if(l[0] == "//"): continue # comment line
    elif(l[0] == "under"): addExp(exps, l[1], int(l[3]), l[4:]) # recordexp
    elif(l[0] == "let"):
     shortcut = l[1]; # ready to enter the braket
     try: shortcut = shortcut[:shortcut.index(':')]
     except ValueError: None;

  return (exps, states);

def preCon(comps, kofe = None, step = 20):
  d = dict();
  for node in comps:
    d[node] = [Bool(node + '_' + str(t)) for t in range(step)];
  return d;

def _get_sublist(l):
  '''
  Generator for non-empty sublists of list l. It will generate sublists in
  non-decreasing size, until it yields the complete list l.
  '''
  for r in range(1, len(l) + 1):
    for combi in combinations(l, r):
      yield list(combi);

def _construct_graph(l, defI = None, allNodes = None):
  '''
  Take a list of directed and signed interactions, return a dictionary of nodes
  with inputs list as value. If defInters are provide, it will be merge in the
  resulting graph.
  '''
  d = dict();
  for i in l:
    appendEntry(d, i[1], (i[0], i[2]));
  if(defI):
    mergeInters(d, defI); # records in defI will be merged in d
  if(allNodes):
    for node in allNodes:
      if(node not in d): d[node] = [];
  return d;

def getGraph(comps, optI, defI):
  '''
  Generator. Yield a non-empty subpgraph from all interactions. Defined
  interactions are bound to be present.
  '''
  l = [];
  for node in optI: # recover dict to list of optional interactions
    l.extend([(i[0], node, i[1]) for i in optI[node]]);
  i = 0; t = 2**len(l);
  complist = comps.keys();
  yield _construct_graph([], defI, complist);
  for sl in _get_sublist(l): # use a subset of optional interactions
      yield _construct_graph(sl, defI, complist);

def sortGraph(graph):
  '''
  '''
  d = dict();
  for node in graph:
    act = filter(lambda x: x[1] > 0, graph[node]);
    rep = filter(lambda x: x[1] < 0, graph[node]);
    d[node] = ([x[0] for x in act], [x[0] for x in rep]);
  return d;

def _compati_func(acrp, funclist):
  (ac, rp) = acrp;
  if(not rp):
    if(not ac): return (-1,) # special function number for isolated node
    elif(len(ac) > 1): return filter(lambda x: x < 2, funclist) or (-1, )
    else: return filter(lambda x: x < 1, funclist) or (-1, )
  if(not ac):
    if(len(rp) > 1): return filter(lambda x: x > 15, funclist) or (-1, )
    else: return filter(lambda x: x > 16, funclist) or (-1,)
  if(len(ac) > 1):
    if(len(rp) > 1): return filter(lambda x: x in (2, 4, 6, 7, 8, 10, 12, 14),
                                   funclist) or (-1,)
    else: return filter(lambda x: x in (2, 4, 6, 7, 10), funclist) or (-1,)
  if(len(rp) > 1):
    return filter(lambda x: x in (2, 4, 14), funclist) or (-1, )
  return filter(lambda x: x in (2, 14), funclist) or (-1, )

def getFunction(comps, sgraph):
  '''
  Generator a combination of all possible combinations of allowed function for
  every node. This generator will also check the compatability of
  rules and input nodes. The yielding is a dictionary, with node:function.
  '''
  nodes = comps.keys();
  rules = dict();
  for n in nodes:
    rules[n] = _compati_func(sgraph[n], comps[n]);
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

def _create_rule(num, act, rep, prec = None, t = None):
  if(num == -1): return False

  if(prec):
    actt = [(prec[node])[t] for node in act]
    rept = [(prec[node])[t] for node in rep]
  else:
    postfix = t and ('_' + str(t)) or '' # if t specified, append _t
    actt = [Bool(a + postfix) for a in act]
    rept = [Bool(r + postfix) for r in rep]

  if(num == 0): return _And(actt)
  elif(num == 1): return _Or(actt)
  elif(num < 4): return And(_Or(actt), Not(_Or(rept))) # at least act, no rep
  elif(num < 6): return And(_And(actt), Not(_And(rept)));
  elif(num == 6): return And(_Or(actt), Not(_And(rept)))
  elif(num == 7): return Or(_And(actt),
                           And(_Or(actt), And(_Or(rept), Not(_And(rept)))))
  elif(num < 10): return _And(actt)
  elif(num < 12): return Or(_And(actt), And(_Or(actt), Not(_Or(rept))))
  elif(num < 14): return Or(_And(actt), And(_Or(actt), Not(_And(rept))))
  elif(num < 16): return _Or(actt)
  elif(num == 16): return And(_Or(rept), Not(_And(rept)))
  elif(num == 17): return Not(_And(rept));

def _with_KOFE(node, kofe):
  if(node in kofe['KO']):
    ko = Bool('KO_' + node);
    if(node in kofe['FE']):
      fe = Bool('FE_' + node);
      return lambda x: Or(fe, And(Not(ko), x)) # both ko and fe
    else: return lambda x: And(Not(ko), x) # only ko
  if(node in kofe['FE']):
    fe = Bool('FE_' + node);
    return lambda x: Or(fe, x)
  else: return lambda x: x; # no kofe

def _add_function(s, fnum, node, kofe, acrp, prec = None, step = 20):
  (ac, rp) = acrp;
  tune = _with_KOFE(node, kofe);
  for t in range(1, step): # 19 transitions from t-1 to t
    s.add(Bool(node + '_' + str(t)) == tune(
               _create_rule(fnum, ac, rp, prec, t = t-1)));

def applyFunctions(s, funcd, kofe, sgraph, prec = None, step = 20):
  '''
  Take the solver, function number, and a sorted graph.
  '''
  for node in funcd:
    _add_function(s, funcd[node], node, kofe, sgraph[node], prec, step);
    
def addConstrains(s, exp, states):
  '''
  exp is a list of tuple, in the form (t, [statenames]). states is the dict
  for shortcut of conditions.
  '''
  for (t, shortcuts) in exp:
    for sc in shortcuts:
      for (node, value) in states[sc]:
        if(node.split('_')[0] in ('KO', 'FE')):
          s.add(Bool(node) == bool(value))
        else: s.add(Bool(node + '_' + str(t)) == bool(value));


### Output Utilities ###
#########################

def printModel(sgraph, funcd):
  s = Solver();
  for node in sgraph:
    (ac, rp) = sgraph[node];
    s.add(Bool(node + "'") == _create_rule(funcd[node], ac, rp));
  print s; print funcd;
