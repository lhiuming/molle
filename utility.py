from z3 import *
from operator import or_
from pprint import pprint

def _sorted_inters(inter_list, sp):
    ''' Sorts the inter_list = [('from', 'to', 'positive'), ...] into a dict,
     where keyvalue is a couple of number tuples, wich integer codes as keys.

    e.g.
    {1: ( (2, 3), (5, 9) )} : species 1 is activated by 2 and 3, and repressed
    by 5 and 9.
    '''
    d = dict([(c, ([], [])) for c in range(len(sp))]) # initialization
    for i in inter_list:
        f, t = i[:2]
        if   'negative' in i: idx = 1
        elif 'positive' in i: idx = 0
        else:
            print 'no +/- assigend to interactions %d'%(inter_list)
            raise(Error)
        tcode, fcode = sp.index(t), sp.index(f)
        d.setdefault(tcode, ([], []))[idx].append(fcode)
    return d

def readModel(f, opmt = True):
    ''' Take a file Object as input, return a tuple of 6 objects:

    species: a tuple of gene name.
    logics : a dict. { gene_name: list_of_allowed_logic_numbers }
    kofe   : a dict. { "FE": list_of_FEable_gene, "KO": list_of_KOable_gene }
    defI   : a dict of defined interations. Processed by _sorted_inters()
    optI   : a dict of optional interactions.
    '''
    species = []
    logics = {}
    kofe = {'KO':[], 'FE':[]}
    def_inters_list = []
    opt_inters_list = []
  
    # read the components line
    for c in f.readline().strip().split(','):
        # get the gene name and +- mark
        if '(' in c : gene_ = c[:c.index('(')].strip()
        else: gene_ = c.strip()
        gene = filter(lambda x: not x in '+-', gene_)
        mark = filter(lambda x: x in '+-', gene_)
        # add to kofe if the gene has mark
        if('+' in mark): kofe['FE'].append(gene)
        if('-' in mark): kofe['KO'].append(gene)
        # record the allowed logics; if no, set to range(18)
        if '(' in c:
            left, right = c.index('('), c.index(')')
            rules = tuple( int(i) for i in c[left+1:right].split() )
        else:
            rules = tuple(range(18))
        logics[gene] = rules
        species.append(gene)

    # read the interaction lines
    total_opt = total_def = 0
    for line in f.readlines():
        l = line.strip().split()
        if(not l): continue # skip empty line
        if 'optional' in l:
            opt_inters_list.append(tuple(l[:3]))
            total_opt += 1
        else:
            def_inters_list.append(tuple(l[:3]))
            total_def += 1
    defI = _sorted_inters(def_inters_list, species)
    optI = _sorted_inters(opt_inters_list, species)

    return (species, logics, kofe, defI, optI)

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

def compati(l, actn, repn):
    ''' Speed up the solving. 
    Not sure with the validicity when actn == 0 of such approach. '''
    if len(l) < 16: return l
    if actn == 0:
        if repn == 0: return (-1, )
        else: # only repressors
            return filter(lambda x: x > 15, l) or (-1, )
    elif repn == 0: # only activator
        return filter(lambda x: x < 2, l) or (-1, )
    else:
        return l

zero = BitVecVal(0, 1)

def Any(bvs):
    return reduce(or_, bvs, zero)

def _concat(bvs):
    if len(bvs) == 1: return bvs[0]
    else: return Concat(bvs)
    
def _create_bit_rule(num, act_list, rep_list, A, R):
    ''' Create the update rule that return bit-vector of length 1. '''
    if num == -1: return BoolVal(False) # special case
    
    # initialization
    if act_list: act = _concat(act_list)
    else: act = A = zero
    if rep_list: rep = _concat(rep_list)
    else: rep = R = zero
    
    # creating result
    if num == 0:
        return And(R == 0, A != 0, A & act == A)
    elif num == 1:
        return And(R == 0, A & act != 0)
       #return And(R == 0, A != 0, A & act != 0)
    elif num == 2:
        return Or( And(R == 0, A != 0, A & act == A),
                   And(R != 0, rep & R == 0, A & act != 0) )
       #return Or( And(R == 0, A != 0, A & act == A),
       #           And(R != 0, A != 0, rep & R == 0, A & act != 0) )
    elif num == 3:
        return And(A & act != 0, rep & R == 0)
    elif num == 4:
        return And( A != 0, A & act == A,
                   Or(R == 0, rep & R != R) )
       #return Or( And(R == 0, A != 0, A & act == A),
       #           And(A != 0, A & act == A, rep & R != R) )
       #return Or( And(R == 0, A != 0, A & act == A),
       #           And(R != 0, A != 0, A & act == A, rep & R != R) )
    elif num == 5:
        return Or( And(R == 0, act & A != 0),
                   And(A != 0, act & A == A, rep & R != R) )
       #return Or( And(R == 0, A != 0, act & A != 0),
       #           And(R != 0, A != 0, act & A == A, rep & R != R) )
    elif num == 6:
        return Or( And(R == 0, A != 0, act & A == A),
                   And(act & A != 0, rep & R != R) )
       #return Or( And(R == 0, A != 0, act & A == A),
       #          And(R != 0, A != 0, act & A != 0, rep & R != R) )
    elif num == 7:
        return Or( And(R == 0, act & A != 0),
                   And(act & A != 0, rep & R != R) )
       #return Or( And(R == 0, A != 0, act & A != 0),
       #           And(R != 0, A != 0, act & A != 0, rep & R != R) )
    elif num == 8:
        return And(A != 0, act & A == A)
       #return Or( And(R == 0, A != 0, act & A == A),
       #           And(R != 0, A != 0, act & A == A) )
    elif num == 9:
        return Or( And(R == 0, act & A != 0),
                   And(R != 0, A != 0, act & A == A) )
       #return Or( And(R == 0, A != 0, act & A != 0),
       #           And(R != 0, A != 0, act & A == A) )
    elif num == 10:
        return Or( And(A != 0, act & A == A),
                   And(R != 0, act & A != 0, rep & R == 0) )
       #return Or( And(R == 0, A != 0, act & A == A),
       #           And(R != 0, A != 0, Or(act & A == A,
       #                                  And(act & A != 0, rep & R == 0))) )
    elif num == 11:
        return Or( And(R == 0, A != 0, act & A != 0),
                   And(R != 0, A != 0, Or(act & A == A,
                                          And(act & A != 0, rep & R == 0))) )
    elif num == 12:
        return Or( And(A != 0, act & A == A),
                   And(act & A != 0, rep & R != R) )
       #return Or( And(R == 0, A != 0, act & A == A),
       #           And(R != 0, A != 0, Or(act & A == A,
       #                                  And(act & A != 0, rep & R != R))) )
    elif num == 13:
        return Or( And(R == 0, A != 0, act & A != 0),
                   And(R != 0, A != 0, Or(act & A == A,
                                          And(act & A != 0, rep & R != R))) )
    elif num == 14:
        return Or( And(R == 0, A != 0, act & A == A),
                   And(R != 0, act & A != 0) )
        #return Or( And(R == 0, A != 0, act & A == A),
        #           And(R != 0, A != 0, act & A != 0) )
    elif num == 15:
        return act & A != 0
        #return Or( And(R == 0, A != 0, act & A != 0),
        #           And(R != 0, A != 0, act & A != 0) )
    elif num == 16:
        return And(A == 0, rep & R != 0, rep & R != R)
        #return And(A == 0, R != 0, rep & R != 0, rep & R != R)
    elif num == 17:
        return And(A == 0, R != 0, rep & R == 0)
    else:
        print "Strange Num"
        raise ValueError

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
    return lambda q, ko, fe: simplify(
        _with_kofe(kofe_index, ko, fe,
                   _create_bit_rule(logic,
                                    [Extract(i,i,q) for i in acts],
                                    [Extract(i,i,q) for i in reps],
                                    A, R)))

def isExpOf2(bvv):
    return len(filter(lambda x: x == '1', bin(bvv.as_long()))) == 1

### Output Utilities ###
#########################
boolf = BoolVal(False)

def conv_time(secs, th = 300):
    if secs > th: return '%.1f min'%( secs / 60 )
    return '%.1f sec'%secs

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
    # simplify is necessary
    return simplify(Extract(i, i, bv)).as_long() == 1

def bv2logic(lbvv, llist):
    ''' convert a bit-vector to a integer, as logic function number.'''
    assert isExpOf2(lbvv)
    lcode = len(bin(lbvv.as_long()).lstrip('0b')) - 1
    return llist[lcode]

def bv2inters(ibvv, ilist, species):
    if simplify(ibvv == 0): return []
    l = ibvv.size() - 1
    return [species[c] for i, c in enumerate(ilist) if checkBit(l-i, ibvv)]

def getDetail(m, A_, R_, L_, species, inters, logics):
    A = {}; R = {}; L = {}
    for c, s in enumerate(species):
        L[s] = bv2logic(m[L_[s]], logics[s])
        if A_[s]: A[s] = bv2inters(m[A_[s]] or zero, inters[c][0], species)
        else: A[s] = []
        if R_[s]: R[s] = bv2inters(m[R_[s]] or zero, inters[c][1], species)
        else: R[s] = []
    return (A, R, L)

def printModel(species, A, R, L, config = True, model = True):
    ''' Print the solved model nicely. '''
    # printing the model
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

from smtplib import SMTP, SMTPAuthenticationError
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

def mailMe(addr, pw, content, title = 'Computation Finished'):
    msg = MIMEMultipart('alternative')
    msg['Subject'] = title
    msg['From'] = msg['To'] = addr
    msg.attach(MIMEText(content, 'plain'))

    server = SMTP('smtp.qq.com')
    try:
        server.login(addr, pw)
        server.sendmail(addr, addr, msg.as_string())
        server.quit()
    except SMTPAuthenticationError:
        print ">> SMTP: login fail with %s:%s"%(addr, pw)
