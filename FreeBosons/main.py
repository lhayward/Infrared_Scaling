import free_boson
import numpy as np
import os.path  #to check if file exists
import sys  #for sys.stdout.flush()
import time

alphaEqTol=1e-10

###############################################################################################
##########################################  arrayEq  ##########################################
###############################################################################################
def arrayEq(a1, a2):
  eq = True
  
  if len(a1) != len(a2):
    eq = False
  else:
    for elem1, elem2 in zip(a1, a2):
      #if elem1 != elem2:
      #Check if the alpha values are equal (within tolerance):
      if abs(elem1-elem2) > alphaEqTol:
        eq = False
        print "%.20f neq %.20f" %(elem1,elem2)
  #end if-else
  
  return eq
####### end arrayEq(a1, a2) function #######

###############################################################################################
########################################  decimalStr  #########################################
###############################################################################################
def decimalStr(num):
  res = str(num)
  length = len(res)
  index = res.find('.')
  if index >= 0:
    res = res[0:index] + '-' + res[(index+1):length]
  return res
####### end decimalStr(num) function #######

###############################################################################################
#########################################  readArray  #########################################
###############################################################################################
def readArray(line):
  start = max( 0, line.find('[') )
  end = min( len(line), line.find(']') )
  return line[start+1:end] 
####### end readArray(line) function #######

###############################################################################################
########################################  readParams  #########################################
###############################################################################################
# Input file should have the following form:
#
# massterm  = ____ (float)
# alpha     = [ ___  ___  ___ ...] (list of floats in square brackets with items separated by spaces)
#
def readParams(filename):
  filename = "input.txt"
  D        = 2
  L_min    = 1
  L_max    = 1
  bc       = 'PBC'
  massterm = 1
  alpha    = [1]
  if os.path.isfile(filename):
    fin = open(filename,'r')
    
    line=fin.readline()
    D = int(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    L_min = int(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    L_max = int(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    bc = line[ max(0,line.find('='))+1:].strip()
    
    line=fin.readline()
    massterm = float(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    alpha = [float(a) for a in readArray(line).split()]
    
    fin.close()
  return D,L_min,L_max,bc,massterm,alpha

###############################################################################################
###########################################  main  ############################################
###############################################################################################

# User settings:
#############################
###Manual input:
# alpha=np.array( np.linspace(0.4,10,49).tolist() + [20,50,100,200,500,1000] )
# #alpha=np.array( [2.2, 2.4] )
# #alpha = np.array( [0.5, 1, 1.5, 2, 2.5, 3])
# massterm = 0.0

###OR: read input from file:
D,L_min,L_max,bc,massterm,alpha = readParams("input.txt")
#############################

BC_options = ['PBC', 'APBC', 'APBCx', 'APBCy', 'OBCDir', 'OBCNeu']
printInDir = False  #whether or not to print into directories names by BC

if bc in BC_options:

  print
  print "dim.  = %d" %D
  print "L_min = %d" %L_min
  print "L_max = %d" %L_max
  print "BC    = " + str(bc)
  print "mass  = %f" %massterm
  print "alpha = " + str(alpha)

  t1 = time.clock()

  total = np.zeros(len(alpha))

  fout_res=[0 for i in alpha]
  for i,n in enumerate(alpha):
    filename = "results_" + bc + "_SVsL_mass" + decimalStr(massterm) + "_alpha" + decimalStr(n) + ".txt"
  
    if printInDir:
      if not(os.path.isdir(bc)):
        os.mkdir(bc)
      filename = bc + '/' + filename
    fout_res[i] = open(filename, 'w')

  for L in range(L_min,L_max+1):
    print
    print "L = %d:" %L
    sys.stdout.flush()

    entropy = free_boson.getEntropy_singleSite(L,bc,alpha,massterm)
    print "  %f" %entropy
    
    # Save result to file
    for i in range(len(alpha)):
      fout_res[i].write("%d %.15f"%(L,entropy[i])+'\n')
      fout_res[i].flush()
  #end loop over orders

  for fout in fout_res:
    fout.close()

  t2 = time.clock()
  print "Total time: " + str(t2-t1) + " sec."
else:
  print "\n*** ERROR: Boundary condition '" + bc + "' not supported ***\n"
