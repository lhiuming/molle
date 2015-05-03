from z3 import *

PREFIX = "examplefiles/";
MODEL = "SimpleFourComponentModel.txt";
EXP = "SingleExperiment.txt";

### Utilities ###
#################
def addInter(d, a, b, p):
  try:
    d[a].append((b, (p == 'positive' and 1 or -1)));
  except KeyError:
    d[a] = [];
    addInter(d, a, b, p);
  return;


### Reading files and Initialization ###
########################################
modelFile = open(PREFIX + MODEL, 'r');

# read components and allowed logic functions
af = modelFile.readline();
comps = dict();
for a in af.strip().split(', '):
  comps[a[0]] = [int(i) for i in a[2:-1].split()];
#print(comps); # debugging

# read interactions and store them as defined or optional interactions
modelFile.readline(); # skip the blank line
defInters = dict();
optInters = dict();
for l in modelFile.readlines(): # loop every entries
  l = l.strip().split();
  if(l[-1] == 'optional'): addInter(optInters, l[0], l[1], l[2]) 
  else: addInter(defInters, l[0], l[1], l[2]); # defined
#print(defInters, optInters); # debugging
modelFile.close();

# read the experiment constrains
expFile = open(PREFIX + EXP, 'r');
expFile.close();



