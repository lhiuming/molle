from z3 import *
from utility import *
from pprint import pprint

### Configurations ###
######################
PREFIX = "examplefiles/"
MODEL = "SimplestModel.txt"
EXP = "constrains.txt"
OUTPUT = "SOLUTIONS.TXT"
STEP = 20
debug = 1 # level 1 of debugging
solutions_limit = 1
interactions_limit = 17


### Reading Files ###
#####################
modelFile = open(PREFIX + MODEL, 'r'); # Read Model File
(comps, kofe, defInters, optInters) = readModel(modelFile);
modelFile.close();

expFile = open(PREFIX + EXP, 'r'); # Read experiment constrains
(exps, states) = readExp(expFile);
expFile.close();

#debuggin
if __name__ == "__main__":
  if(debug == 2):
    print("The components are: \n");
    pprint(comps);
    print("\nThe defined and optional interactions: \n");
    pprint(defInters);
    pprint(optInters);
    print("\nThe Experiment Constrains: \n");
    pprint(exps);
    pprint(states + '\n\n');


### Modeling ###
################

count = 0; # count model numbers
solution_count = 0; # count solutions
precon = preCon(comps, kofe, step = STEP); # preconstruct z3 Bool variable

for graph in getGraph(comps, optInters, defInters, interactions_limit):
  sgraph = sortGraph(graph); # sort the graph dict into a more compact form
  if(debug):
    print ">>> The current graph: "
    pprint(sgraph)

  for funcDict in getFunction(comps, sgraph):
    count += 1; # couting solutions
    if(count % 100 == 0):
      print ">>> Verifying the %d model:" %count;
      if(debug): printModel(sgraph, funcDict);

    s = Solver();
    applyFunctions(s, funcDict, kofe, sgraph, precon, STEP); # construct model
    for name in exps: # test experiment constrains seperately
      s.push();
      addConstrains(s, exps[name], states);
      if(s.check() == unsat): break; # must satisfy all experiment constrains
      s.pop();

    if(s.check() == sat): # all constrains satisfied
      solution_count += 1;
      outputModel(sgraph, funcDict, OUTPUT, solution_count);
      break; # one combination of funcionts is enough for one graph

  if(solution_count >= solutions_limit): break;

if(solution_count == 0)
  print ">>> No solutions possible, after verifying %d models.\n" %count
