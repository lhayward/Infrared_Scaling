import free_fermion_2D
import numpy as np
import argparse

parser=argparse.ArgumentParser(description="Full diagonalization of free fermion Hamiltonian in 2D for squares")

parser.add_argument('--L','-l',type=int,nargs='+',help='System size, range possible',required=True)
parser.add_argument('--Ratio','-r',type=int,default=2,help='Ratio between linear extent of subsystem A and total system (--L), default=2')
parser.add_argument('--BC','-b',default='obc',help='Boundary conditions (pbc, obc, apc), default=obc')
parser.add_argument('--Mass','-m',type=float,default=1.0,help='Mass term, default=1.0')

alpha=[1.0]
args=parser.parse_args()

if len(args.L)==1:
    L=[args.L]
else:
    L=range(args.L[0],args.L[1]+1,args.L[2])

ratio=args.Ratio

print "L_min = %d" %min(L)
print "L_max = %d" %max(L)
print "BC   = " + str(args.BC)
print "mass  = %f" %args.Mass
print "alpha = " + str(alpha)

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

def RegionA(L,LA,offX=0,offY=0,obc=True):
    A=np.zeros((2*L,L),dtype=bool)
    #A[L/2-LA/2+offX:L/2-LA/2+LA+offX,2*(L/2-LA/2+offY):2*(L/2-LA/2+LA+offY)]=True
    #print A.T
    
    if obc==True:
        A[2*(L/2-LA/2+offX):2*(L/2-LA/2+LA+offX),L/2-LA/2+offY:L/2-LA/2+LA+offY]=True
    else:
        A[0:2*LA,0:LA]=True
    return A

def getHamiltonian(Lx,Ly,bc='pbc',massterm=1.0): #accepts also 'obc'
    Ns=Lx*Ly
    local_diag=np.diag([massterm,-massterm]) #True value on the diagonal is 2*m, but we achieve this by adding the hermitian conjugate.
    hopping_x_p=np.matrix(([-0.5*massterm, -0.5*1j],[-0.5*1j, 0.5*massterm]), dtype='complex')
    hopping_y_p=np.array(([-0.5*massterm, -0.5],[0.5, 0.5*massterm]))
    
    diag=np.kron(np.eye(Ns),local_diag)
    if bc=='pbc':
        diag_x=np.roll(np.eye(Lx),1,axis=1)
        hopping_x=np.kron(diag_x,hopping_x_p)
        offdiag_x=np.kron(np.eye(Ly),hopping_x)
    
        offdiag_y=np.kron(np.roll(np.eye(Ns),Lx,axis=1),hopping_y_p)
    
    if bc=='yapbc':
        diag_x=np.roll(np.eye(Lx),1,axis=1)
        hopping_x=np.kron(diag_x,hopping_x_p)
        offdiag_x=np.kron(np.eye(Ly),hopping_x)
    
        offdiag_y=np.kron(np.diag([1.]*(Ns-Lx),Lx),hopping_y_p)
        offdiag_y+=np.kron(np.diag([1.]*Lx,(Ns-Lx)),-hopping_y_p.T)
    
    elif bc=='obc':
        diag_x=np.diag([1.]*(Lx-1),1)
        hopping_x=np.kron(diag_x,hopping_x_p)
        offdiag_x=np.kron(np.eye(Ly),hopping_x)
        
        offdiag_y=np.kron(np.diag([1.]*(Ns-Lx),Lx),hopping_y_p)
    
    upper_tri=diag+offdiag_x+offdiag_y
    Ham=upper_tri+upper_tri.getH()
    
    
    return Ham


np.set_printoptions(linewidth=200)

if args.BC == 'obc':
    filename = "results_OBC_SVsL_mass" + decimalStr(args.Mass) + "_alpha1-0.txt"
    obc = True
else:
    filename = "results_PBC_SVsL_mass" + decimalStr(args.Mass) + "_alpha1-0.txt"
    obc = False
fout = open(filename, 'w')

for l in L:
    print "L = %d:" %l
    entropy=free_fermion_2D.getEntropy((l,l),alpha,RegionA(l,1,0,0,obc),getHamiltonian(l,l,args.BC,massterm=args.Mass))[0]
    #print l/ratio, free_fermion_2D.getEntropy((l,l),alpha,RegionA(l,l/ratio,False),getHamiltonian(l,l,args.BC,massterm=args.Mass))[0]
    #print l, free_fermion_2D.getEntropy((l,l),[1.0],RegionA(l,1),getHamiltonian(l,l,args.BC))[0]
    print "  %f" %entropy
    fout.write("%d %.15f"%(l,entropy)+'\n')

fout.close()
