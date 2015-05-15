### Configurations ###
######################

# Running Setting
processN = 4 # processes number
debug = 1 # debugging level

# computation settings
solutions_limit = 1 # how many solution you wish to find
interactions_limit = 17 # set to 17 for the minimal model
graphs_limit = 1 # limits graph numbers. Don't used normally

# Input and output files
OUTPUT = "solutions.out" # file name for solution output
PREFIX = "examplefiles/"
INPUT = { 'ABCD_test': ( "SimpleFourComponentModel.txt",
                         "CertainInteractionRequired.txt" ),
          'minimal_test': ( "custom.txt",
                            "UltimateConstrains.txt" ),
          "find_minimal_combination": ( "SimplestModel.txt",
                                        "UltimateConstrains.txt" ),
          "find_minimal_model":
              ( "PearsonThreshold792WithOct4Sox2Interaction.txt",
                "UltimateConstrains.txt" )
        }
MODEL, EXP = INPUT['minimal_test']

# Model configuration
STEP = 20 # trajactory length
use_compact = True # use compact list of allowed function


### Modelling #########################################
#######################################################
from z3 import *

from multiprocessing import Process, Pipe
from Queue import Queue

from utility import *
from pprint import pprint
from time import sleep, strftime

def now(): return strftime("%d %b %H:%M >>")

def checkGraph(comps, kofe, graph, exps, states, precon, compact, gnum, conn):
  '''
  Worker function for children processes
  '''
  sgraph = sortGraph(graph); # sort the graph dict into a more compact form  
  funcds = getFunction(comps, sgraph, compact)

  # reporting
  ftotal, rules = funcds.next() # the first yield is configuration information
  if(debug):
    print now(), "Start: Graph %d have %d combinations." %(gnum, ftotal)
    if(debug > 1):
      print "> The graph: ", sgraph
      print "> The funcion lists: ", rules

  # check all combinations
  model_count = 0
  for funcd in funcds:
    # reporting
    model_count += 1
    if(debug and model_count % 500  == 0):
      print now(), "Graph %d has verified %d/%d ." %(gnum, model_count, ftotal)
    # buld up the model
    s = Solver()
    applyFunctions(s, funcd, kofe, sgraph, precon, STEP)
    # chekc against experimental constrains
    passed = True
    for name in exps: # test against all constrains
      s.push()
      addConstrains(s, exps[name], states, precon)
      if(s.check() == unsat):
        passed = False
        break
      s.pop()
    # if get a passed model, send this solution
    if(passed):
      conn.send( ( sgraph, funcd) )
      break;    
  # if no solution found after checking all combinations
  if(not passed): conn.send(0)

  # send the graph number, and close the connection
  conn.send((gnum, model_count))
  conn.close()

def collectResult(p, conn, solutions_count, model_num):
  '''
  Collect result from process p, and output the result if it got a solution.
  '''
  result = conn.recv()
  gnum, mcount = conn.recv()
  model_num += mcount
  p.join() # end the children process
  if(result): # if a solution is received
    solutions_count += 1
    print now(),"The %dth solution found at Graph %d." %(solutions_count, gnum)
    outputModel(*result, fileName = OUTPUT, count = solutions_count)
  print now(),"Verified %d models in total" %model_num
  return solutions_count, model_num

### Main Process
if __name__ == "__main__":
  # Read model files and experiment constrains file
  modelFile = open(PREFIX + MODEL, 'r')
  (comps, kofe, defInters, optInters) = readModel(modelFile); modelFile.close()
  expFile = open(PREFIX + EXP, 'r')
  (exps, states) = readExp(expFile); expFile.close();
  
  # Pre-construct z3 Bool objects. May save computation
  precon = preCon(comps, kofe, step = STEP)
  
  # multiprocess computation
  solutions_count = graphs_count = model_count = 0 # intialize the counters
  workers = Queue(processN) # a queue for workers processes
  for graph in getGraph(comps, optInters, defInters, interactions_limit):
    graphs_count += 1
    
    # make sure there is avaliable site in the queue for a more worker
    while(workers.full()): # wait for a empty site
      p, conn = workers.get()
      if(conn.poll(5)): # get finishing signal, waiting for 5s
        solutions_count, model_count = collectResult(p, conn, solutions_count,
                                                     model_count)
        break
      else: workers.put( (p, conn) ) # not finished yet, put it back

    if(solutions_limit and solutions_count >= solutions_limit):
      print now(),"solutions limit reached. No more graph"
      break
      
    # create a new worker for the new graph
    parent_conn, child_conn = Pipe() # pipe for comunication between processes
    p = Process(target=checkGraph,
                args=(comps, kofe, graph, exps, states, precon, use_compact,
                      graphs_count, child_conn))
    p.start()
    workers.put( (p, parent_conn) ) # put the new worker in queue
    
    if(graphs_limit and graphs_count >= graphs_limit):
      print now(),"graphs limit reached. No more graph"
      break
    
  # collect the rest works
  while(not workers.empty()):
    (p, conn) = workers.get()
    if(solutions_limit and solutions_count >= solutions_limit):
      p.terminate(); conn.close()      
    elif(conn.poll(5)):
      solutions_count, model_count = collectResult(p, conn, solutions_count,
                                                   model_count)
    else: workers.put( (p, conn) )

  # final reporting
  if(solutions_count == 0):
    print now(),"No solutions found, after verifying %d models on %d graphs.\n"\
      %(model_count, graphs_count)
  else:
    print now(),"Done. %d solutions found, after verifying %d models"\
          " on %d graphs.\n" %(solutions_count, model_count, graphs_count)
