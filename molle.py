from z3 import *
from utility import *

from pprint import pprint

### Configurations ###
######################
PREFIX = "examplefiles/";
MODELFILES = ("SimpleFourComponentModel.txt",
              "PearsonThreshold792.txt");
EXPFILES = ("CertainInteractionRequired.txt", "NoSolutionsPossible.txt",
            "ExperimentalConstraints.txt");
MODEL = MODELFILES[1];
EXP = EXPFILES[2];
STEP = 20;
debug = True;
limit = 10;


### Reading Files ###
#####################

# Read Model File
modelFile = open(PREFIX + MODEL, 'r');
(comps, kofe, defInters, optInters) = readModel(modelFile);
modelFile.close();

# Read experiment constrains
expFile = open(PREFIX + EXP, 'r');
(exps, states) = readExp(expFile);
expFile.close();


### Modeling ###
################

if __name__ == "__main__":
  if(debug):
    print("The components are: ");
    pprint(comps);
    print("\nThe defined and optional interactions: ");
    pprint(defInters);
    pprint(optInters);
    print("\nThe Experiment Constrains: ");
    pprint(exps);
    pprint(states);

  count = 0;
  count_solution = 0;
  graphs = getGraph(comps, optInters, defInters);
  functions = getFunction(comps);

  for graph in graphs:
    sgraph = sortGraph(graph);

    for funcDict in functions:
      count += 1;
      if(count % 100 == 0): print "doing the %d model..." %count;

      s = Solver();
      applyFunctions(s, funcDict, kofe, sgraph, STEP);
      for name in exps:
        s.push();
        addConstrains(s, exps[name], states);
        if(s.check() == unsat): break;
        s.pop();

      if(s.check() == sat):
        count_solution += 1;
        printModel(sgraph, funcDict);

      if(count_solution >= limit): break;
    if(count_solution >= limit):break;

  if(count_solution == 0):
    print "No solutions possible, after %d times search." %count;    

