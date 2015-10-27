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
  L_min    = 1
  L_max    = 1
  per      = False
  massterm = 1
  alpha    = [1]
  if os.path.isfile(filename):
    fin = open(filename,'r')
    
    line=fin.readline()
    L_min = int(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    L_max = int(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    per = bool(int(line[ max(0,line.find('='))+1:]))
    
    line=fin.readline()
    massterm = float(line[ max(0,line.find('='))+1:])
    
    line=fin.readline()
    alpha = [float(a) for a in readArray(line).split()]
    
    fin.close()
  return L_min,L_max,per,massterm,alpha

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
L_min,L_max,per,massterm,alpha = readParams("input.txt")
#############################

print
print "L_min = %d" %L_min
print "L_max = %d" %L_max
print "PBC   = " + str(per)
print "mass  = %f" %massterm
print "alpha = " + str(alpha)

t1 = time.clock()

total = np.zeros(len(alpha))

fout_res=[0 for i in alpha]
for i,n in enumerate(alpha):
  filename = "results_SVsL_mass" + decimalStr(massterm) + "_alpha" + decimalStr(n) + ".txt"
  fout_res[i] = open(filename, 'w')

for L in range(L_min,L_max+1):
  print
  print "L = %d:" %L
#   for Lx,Ly,Lz in order.clusters(ord):
#     curr_clust_name = clust_weight.clust_name(Lx,Ly,Lz)
#     print "  " + curr_clust_name
#     sys.stdout.flush()
#     
#     w = clust_weight.weight(Lx,Ly,Lz,w,alpha,massterm) # performs cluster weight calculations
#     #Embedding factor:
#     ef = 6
#     if Lx==Ly and Ly==Lz: ef = 1
#     elif Lx==Ly or Lx==Lz or Ly==Lz: ef = 3
#     
#     total = total + ef*w[curr_clust_name]
#   #end loop over loop over clusters
#     
#   # Save result to file
#   for i in range(len(alpha)):
#     fout_res[i].write("%d %.15f"%(ord,total[i])+'\n')
#     fout_res[i].flush()
#end loop over orders

for fout in fout_res:
  fout.close()

t2 = time.clock()
print "Total time: " + str(t2-t1) + " sec."
