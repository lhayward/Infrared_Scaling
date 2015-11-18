import numpy as np
import pprint

###############################################################################################
########################################  getEntropy  #########################################
###############################################################################################
def getEntropy((Lx,Ly),alpha,regA,Hamiltonian):
    sitesA = []
    for x in range(2*Lx):
        for y in range(Ly):
            if regA[x,y] == True: sitesA.append(site((2*Lx,Ly),x,y))
    NsA=len(sitesA)
    #Probably np.flatten does it
    sitesA.sort()
    
    Eval,Evec = np.linalg.eigh(Hamiltonian) 
    U=np.matrix(Evec)
    #print Eval
    
    '''
    GSOccup=np.kron(([1,0],[0,0]),np.eye(Lx*Ly))
    C= U*GSOccup*U.getH()
    '''
    GSU=U[:,0:Lx*Ly]
    C=GSU*GSU.getH()
    
    Cred=C[np.ix_(sitesA,sitesA)]
    
    Ev = abs(np.linalg.eigvals(Cred).real)
    #print Ev
    
    Sn = np.zeros(len(alpha))
    #Sn = np.ones(len(alpha))
    
    vonNeumannS=0.0
    vonNeumannS = - sum(Ev*np.log(Ev+1e-15)) - sum((1.-Ev)*np.log(abs(1.-Ev+1e-15)))
    for i,n in enumerate(alpha):
        if n==1:
            Sn[i]=vonNeumannS
        else:
            Sn[i]=(1./(1.-n))*sum(np.log((Ev+1e-15)**n + (abs(1-Ev+1e-15))**n))
            #Sn[i]=(1./(1.-n))*(Sn[i])
    #print Sn
    return Sn
#........................................END getEntropy........................................

###############################################################################################
########################################  getRegionA  #########################################
###############################################################################################

def getRegionA((Lx,Ly),(x0,y0),bipart): #all factors of two due to fermion spin 
    regA = np.zeros( (2*Lx,Ly), dtype='bool' ) #doubling in x direction only!
    if bipart==1:
        regA[:2*x0,:]=True
    elif bipart==2:
        regA[:,:y0]=True
    elif bipart==3:
        regA[2*x0:,y0:]=True
    else:
        regA[:2*x0,:y0]=True

    if( (regA==True).sum() > (2*Lx*Ly/2) ): 
        regA = np.logical_not(regA)
    
    print regA

    return regA

###############################################################################################
###########################################  site  ############################################
###############################################################################################
def site((Lx,Ly),x,y):
  return x + (y*Lx)  #convert (x,y) pair to a single site number
