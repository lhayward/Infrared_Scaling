import numpy as np
#from scipy import linalg
#import time
#start_time = time.time()

#Ev_tol = 10e-16

###############################################################################################
################################  getEntropy_singleSite_2D_K  #################################
###############################################################################################
def getEntropy_singleSite_2D_K(L,bc,alpha,massterm):
  #PBC or not along the x and y directions:
  perx = pery = False
  if bc == 'PBC':
    perx = pery = True
  Lx = Ly = L
  
  Ns   = Lx * Ly
  Ltup = (Lx,Ly)
  
  K = np.zeros((Ns, Ns))
  #Loop to calculate the matrix K:
  for x in range(Lx):
    for y in range(Ly):
      K[site(Ltup,x,y),site(Ltup,x,y)]= 4.0 + (massterm ** (2))
      xp = (x+1)%Lx
      yp = (y+1)%Ly
      
      #Adjust onsite term for Neumann BC:
      if bc == 'OBCNeu':
        if xp < x or x==0:  #if we are on the left or right boundary
          K[site(Ltup,x,y),site(Ltup,x,y)] = K[site(Ltup,x,y),site(Ltup,x,y)] - 1
        if yp < y or y==0:  #if we are on the bottom or top boundary
          K[site(Ltup,x,y),site(Ltup,x,y)] = K[site(Ltup,x,y),site(Ltup,x,y)] - 1
      
      if (xp > x) or perx:
        K[site(Ltup,xp,y), site(Ltup,x, y)] = K[site(Ltup,xp,y), site(Ltup,x, y)] - 1.0 
        K[site(Ltup,x, y), site(Ltup,xp,y)] = K[site(Ltup,x, y), site(Ltup,xp,y)] - 1.0
      if (yp > y) or pery:
        K[site(Ltup,x,yp), site(Ltup,x,y )] = K[site(Ltup,x,yp), site(Ltup,x,y )] - 1.0 
        K[site(Ltup,x,y ), site(Ltup,x,yp)] = K[site(Ltup,x,y ), site(Ltup,x,yp)] - 1.0
  
  #use eigh because we know the matrix is symmetric
  Eval,Evec = np.linalg.eigh(K) #, b=None, left=False, right=True, overwrite_a=False, overwrite_b=False, check_finite=True)
  #print Eval
  #print Evec
  P = 1./2. * np.matrix(Evec) * np.matrix(np.diag(np.sqrt(abs(Eval)))) * np.matrix(Evec.T)
  X = 1./2. * np.matrix(Evec) * np.matrix(np.diag(1. / np.sqrt(abs(Eval)))) * np.matrix(Evec.T)
  
  ###
  #sitesA = [i for i in range(len(regA)) if regA[i]]
  sitesA = [ site(Ltup, int(Lx)/2, int(Ly)/2) ]
  #sitesA = [ site(Ltup, 0, 0) ]
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
#................................END getEntropy_singleSite_2D_K................................

###############################################################################################
###################################  getXAPA_singleSite_2D  ###################################
###############################################################################################
def getXAPA_singleSite_2D(L,bc,alpha,massterm):
  Lx = Ly = L
  x_A = int(Lx)/2
  y_A = int(Ly)/2
  
  def omega_2D(qx,qy):
    return np.sqrt( 4*np.sin(qx/2.0)**2 + 4*np.sin(qy/2.0)**2 + massterm**2 )
  
  XA = 0
  PA = 0
  print bc
  if bc in ['PBC', 'APBC', 'APBCy', 'APBCx']:   #PBC or APBC
    for nx in range(Lx):
      for ny in range(Ly):
        if bc == 'PBC':
          kx = 2*nx*np.pi/(Lx)
          ky = 2*ny*np.pi/(Ly)
        elif bc == 'APBC':
          kx = (2*nx + 1)*np.pi/(Lx)
          ky = (2*ny + 1)*np.pi/(Ly)
        elif bc == 'APBCx':
          kx = (2*nx + 1)*np.pi/(Lx)
          ky = 2*ny*np.pi/(Ly)
        else: #APBCy (should be the same as APBCx when Lx=Ly for single site)
          kx = 2*nx*np.pi/(Lx)
          ky = (2*ny + 1)*np.pi/(Ly)
          
        #Set x_A = y_A = 0:
        XA = XA + (1.0/omega_2D(kx,ky))
        PA = PA + omega_2D(kx,ky)
      #end loop over ny
    #end loop over nx
      
    XA = 1.0/(2.0*Lx*Ly)*XA
    PA = 1.0/(2.0*Lx*Ly)*PA
  elif bc == 'OBCDir':   #Dirichlet OBC
    for nx in range(1,Lx+1):
      for ny in range(1,Ly+1):
        qx = nx*np.pi/(Lx+1)
        qy = ny*np.pi/(Ly+1)
        #use (x_A + 1) because Sharmistha's paper uses indices starting at 1
        XA = XA + np.sin(qx*(x_A+1))**2*np.sin(qy*(y_A+1))**2/omega(qx,qy)
        PA = PA + np.sin(qx*(x_A+1))**2*np.sin(qy*(y_A+1))**2*omega(qx,qy)
    XA = 2.0/((Lx+1)*(Ly+1))*XA
    PA = 2.0/((Lx+1)*(Ly+1))*PA
  #end OBC case
  else:
    print "\n*** ERROR: In 2D, the boundary condition '" + bc + "' is not accepted.***\n"
  
  return XA,PA

###############################################################################################
###################################  getEntropy_singleSite  ###################################
###############################################################################################
# This method uses the known formulas for the correlation functions (for periodic or Dirichlet
# boundary conditions). 
def getEntropy_singleSite(L,bc,alpha,massterm):
  XA, PA = getXAPA_singleSite_2D(L,bc,alpha,massterm)
  CA = np.sqrt(XA*PA)
  
  Sn = np.zeros(len(alpha))
  for i,n in enumerate(alpha):
    if CA > 1.0/2:  #numerical issues when CA=1/2 exactly, but Sn should be zero in this case
      if n == 1:
        Sn[i] = (CA+1./2)*np.log(CA+1./2.) - (CA-1./2.)*np.log(CA-1./2)
      else:
        Sn[i] = 1.0/(n-1.0)*np.sum( np.log( (CA+1./2)**n - (CA-1./2.)**n ) )
    else:
      Sn[i] = 0
  
  return Sn
#................................END getEntropy_singleSite_2D..................................



###############################################################################################
###########################################  site  ############################################
###############################################################################################
def site((Lx,Ly),x,y):
  return x + (y*Lx) #convert (x,y) pair to a single site number
