from z3 import *
from utility import *
from pprint import pprint

### Configurations ###
######################
PREFIX = "examplefiles/"
MODEL = "SimplestModel.txt"
EXP = "constrains.txt"
#MODEL = "SimpleFourComponentModel.txt"
#EXP = "CertainInteractionRequired.txt"
OUTPUT = "SOLUTIONS.TXT"
STEP = 20
debug = 1
solutions_limit = 1
interactions_limit = 17
graph_stop = 1


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
    print("The components are:");
    pprint(comps);
    print("\nThe defined and optional interactions:");
    pprint(defInters);
    pprint(optInters);
    print("\nThe Experiment Constrains:");
    pprint(exps);
    pprint(states);


### Modeling ###
################

count = 0; # count model numbers
solution_count = 0; # count solutions
graph_count = 0;
precon = preCon(comps, kofe, step = STEP); # preconstruct z3 Bool variable
s = Solver()

for graph in getGraph(comps, optInters, defInters, interactions_limit):
  sgraph = sortGraph(graph); # sort the graph dict into a more compact form
  graph_count += 1
  if(debug):
    print ">>> The current graph: "
    pprint(sgraph)

  for funcDict in getFunction(comps, sgraph):
    count += 1; # couting solutions
    if(count % 100 == 0):
      print ">>> Verifying the %d model:" %count;
      if(debug): printModel(sgraph, funcDict);

    s.push()
    applyFunctions(s, funcDict, kofe, sgraph, precon, STEP); # construct model
    for name in exps: # test experiment constrains seperately
      s.push();
      addConstrains(s, exps[name], states);
      if(s.check() == unsat): break; # must satisfy all experiment constrains
      s.pop();

    if(s.check() == sat): # all constrains satisfied
      solution_count += 1;
      outputModel(sgraph, funcDict, OUTPUT, solution_count);
      s.pop()
      break; # one combination of funcionts is enough for one graph
    s.pop()

  if(solution_count >= solutions_limit or graph_count >= graph_stop): break;

if(solution_count == 0):
  print ">>> No solutions possible, after verifying %d models.\n" %count
else:
  print ">>> Done. %d solutions were found, after verifying %d models.\n" \
        %(solution_count, count)
