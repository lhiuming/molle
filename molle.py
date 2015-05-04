from z3 import *
from utility import *
from pprint import pprint

### Configurations ###
######################
PREFIX = "examplefiles/";
MODEL = "SimpleFourComponentModel.txt";
EXP = "SingleExperiment.txt";
debug = True;


### Reading Files ###
#####################

# Read Model File
modelFile = open(PREFIX + MODEL, 'r');
(comps, defInters, optInters) = readModel(modelFile);
modelFile.close();

# Read experiment constrains
expFile = open(PREFIX + EXP, 'r');
(exps, states) = readExp(expFile);
expFile.close();



### Modeling ###
################




### Debugging ###
#################

if(debug):
  print("The components are: ");
  pprint(comps);
  print("\nThe defined and optional interactions: ");
  pprint(defInters);
  pprint(optInters);

  print("\nThe Experiment Constrains: ");
  pprint(exps);
  pprint(states);
