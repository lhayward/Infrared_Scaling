import numpy as np
#from scipy import linalg
#import time
#start_time = time.time()

#Ev_tol = 10e-16

###############################################################################################
########################################  getEntropy  #########################################
###############################################################################################
def getEntropy(L,per,alpha,massterm):
  perx = pery = per  #PBC or not along the x and y directions
  Lx = Ly = L
  
  Ns   = Lx * Ly
  Ltup = (Lx,Ly,Lz)
  
  K = np.zeros((Ns, Ns))
  #Loop to calculate the matrix K:
  for x in range(Lx):
    for y in range(Ly):
      K[site(Ltup,x,y),site(Ltup,x,y)]= 4.0 + (massterm ** (2))
      xp = (x+1)%Lx
      yp = (y+1)%Ly
      if (xp > x) or perx:
        K[site(Ltup,xp,y), site(Ltup,x, y)] = -1.0 
        K[site(Ltup,x, y), site(Ltup,xp,y)] = -1.0
      if (yp > y) or pery:
        K[site(Ltup,x,yp), site(Ltup,x,y )] = -1.0 
        K[site(Ltup,x,y ), site(Ltup,x,yp)] = -1.0
  
  #use eigh because we know the matrix is symmetric
  Eval,Evec = np.linalg.eigh(K) #, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True)
  #print Evec
  P = 1./2. * np.matrix(Evec) * np.matrix(np.diag(np.sqrt(Eval))) * np.matrix(Evec.T)
  X = 1./2. * np.matrix(Evec) * np.matrix(np.diag(1. / np.sqrt(Eval))) * np.matrix(Evec.T)
  
  ###
  sitesA = [i for i in range(len(regA)) if regA[i]]   
  Pred = P[sitesA][:,sitesA]
  Xred = X[sitesA][:,sitesA]
  
  Csquared = np.matrix(Xred)*np.matrix(Pred)
  
  Ev = np.sqrt(np.linalg.eigvals(Csquared))
  
  Sn = np.zeros(len(alpha))
  Ev_new = np.array([e.real for e in Ev if (e.real - 1./2.)>0])
  for i,n in enumerate(alpha):
    if n == 1:
      Sn[i] = np.sum( (Ev_new+1./2)*np.log(Ev_new+1./2.) - (Ev_new-1./2.)*np.log(Ev_new-1./2) )
    else:
      Sn[i] = 1.0/(n-1.0)*np.sum( np.log( (Ev_new+1./2)**n - (Ev_new-1./2.)**n ) )
  
  return Sn
  
#   #result=0
#   result = np.zeros(len(alpha))
#   #loop over all possible locations of the corner:
#   for x0 in range(1,Lx):
#     for y0 in range(1,Ly):
#       for z0 in range(1,Lz):
#         r0 = (x0,y0,z0)
#         
#         ################################# CORNER CONTRIBUTION #################################
#         #corners:
#         result += 1.0/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.geq,cmp.geq),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt, cmp.lt, cmp.geq),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.lt, cmp.lt ),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt, cmp.geq,cmp.lt ),X,P))
#         #edges:
#         result -= 1.0/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.geq,cmp.any),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt, cmp.lt ,cmp.any),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.any,cmp.geq),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.lt ,cmp.any,cmp.lt ),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.geq,cmp.geq),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.lt ,cmp.lt ),X,P))
#         #planes:
#         result += 1.0/4.0*(getEntropy(L,alpha,getRegionA(L,r0,cmp.geq,cmp.any,cmp.any),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.geq,cmp.any),X,P)
#                          + getEntropy(L,alpha,getRegionA(L,r0,cmp.any,cmp.any,cmp.geq),X,P))
#                          
#       #loop over z0
#     #loop over y0
#   #loop over x0
#   return result
#.......................................END getCornerEnt.......................................
  

# ###############################################################################################
# ########################################  getEntropy  #########################################
# ###############################################################################################
# def getEntropy((Lx,Ly,Lz),alpha,regA,X,P):
#   sitesA = [i for i in range(len(regA)) if regA[i]]   
#   Pred = P[sitesA][:,sitesA]
#   Xred = X[sitesA][:,sitesA]
#   
#   Csquared = np.matrix(Xred)*np.matrix(Pred)
#   
#   Ev = np.sqrt(np.linalg.eigvals(Csquared))
#   
#   Sn = np.zeros(len(alpha))
#   Ev_new = np.array([e.real for e in Ev if (e.real - 1./2.)>0])
#   for i,n in enumerate(alpha):
#     if n == 1:
#       Sn[i] = np.sum( (Ev_new+1./2)*np.log(Ev_new+1./2.) - (Ev_new-1./2.)*np.log(Ev_new-1./2) )
#     else:
#       Sn[i] = 1.0/(n-1.0)*np.sum( np.log( (Ev_new+1./2)**n - (Ev_new-1./2.)**n ) )
#   
#   return Sn
# #........................................END getEntropy........................................

###############################################################################################
########################################  getRegionA  #########################################
###############################################################################################
def getRegionA((Lx,Ly,Lz),(x0,y0,z0),fx,fy,fz):
  NTot=Lx*Ly*Lz
  regA = np.zeros( NTot, dtype='bool' )
  
  for x in range(Lx):
    for y in range(Ly):
      for z in range(Lz):
        if (fx(x,x0) and fy(y,y0) and fz(z,z0)):
          regA[site((Lx,Ly,Lz),x,y,z)] = True
        #endif
      #end for loop over z
    #end for loop over y
  # end for loop over x
  
  #if the size of region A is more than half of the lattice, swap regions A and B:
  if( (regA==True).sum() > (NTot/2.0) ): regA = np.logical_not(regA)
  
  return regA

###############################################################################################
###########################################  site  ############################################
###############################################################################################
def site((Lx,Ly),x,y):
  return x + (y*Lx) #convert (x,y) pair to a single site number
