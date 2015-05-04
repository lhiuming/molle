from z3 import *

def __init__():
  return;

def appendEntry(d, k, entry):
  try: d[k].append(entry);
  except KeyError:
    d[k] = [];
    appendEntry(d, k, entry);

def addInter(d, a, b, p):
  appendEntry(d, b, (a,(p == 'positive') and 1 or -1));

def readModel(f):
  '''
  Take a file as input, return three dicts, which are allowed rules for 
  components,  defined interactions and optional interactions, seperately.
  The entry of interactions dict collect "from" nodes.
  '''
  comps = dict();
  defInters = dict();
  optInters = dict();

  for c in f.readline().strip().split(', '): # read the components line
    comps[c[0]] = [int(i) for i in c[2:-1].split()]; # store allowed rules num
  f.readline(); # skip the blank line
  for l in f.readlines(): # loop every line for interaction
    l = l.strip().split();
    if(l[-1] == 'optional'): addInter(optInters, l[0], l[1], l[2])
    else: addInter(defInters, l[0], l[1], l[2]); # defined interaction

  return (comps, defInters, optInters);

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
    l = l.split();
    if(not l): continue; # skip empty line
    if(shortcut): # inside the braket { }
      if(l[0] == '{'): continue # skip left bracket
      elif(l[0] == '};'): shortcut = '' # exit the braket;
      else:
        addState(states, shortcut, l[0], int(l[2])); # record configuration
    elif(l[0] == "//"): continue # comment line
    elif(l[0] == "under"): addExp(exps, l[1], int(l[3]), l[4:]) # record states
    elif(l[0] == "let"): shortcut = l[1]; # ready to enter the braket

  return (exps, states);
