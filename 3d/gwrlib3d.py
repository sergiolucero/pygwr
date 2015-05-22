import os
from time import time
from numpy import *
from scipy import sparse, linalg, io
from scipy.sparse import spdiags, csc_matrix, coo_matrix
from scipy.sparse.linalg import linsolve
from matplotlib import pyplot, colors, ticker, use
set_printoptions(suppress=True,precision=4)
from scipy.optimize import fmin_slsqp
import matplotlib
#########################################################################
##############         3D code starts here!               ###############
#########################################################################
class Hydrology3d(object):
    """ a hydrology class """
    def __init__ (self,n0scen,Hwb,Dx,Dy,Dz,WL):
        self.n0scen = n0scen; self.Hwb=Hwb;
        self.Dx=Dx; self.Dy=Dy; self.Dz=Dz;self.WL=WL
#########################################################################
def startph3d(nvec,nscen,nwells):
    t0=time()
    nx,ny,nz=nvec
    print '[PH-3D Started][nxyz=(%d,%d,%d) nscen=%d]' %(nx,ny,nz,nscen)
    LB = -1; UB = 1; rho = 1;
    W0=zeros((nwells,1));   qhat=zeros((nwells,1))
    save('work/refcard3d',[rho,nscen,nwells,nx,ny,nz,LB,UB])
    cmd = 'mpiexec -n %d python presolve3d.py' %(nscen)
    os.system(cmd)
    for scen in range(nscen):
        fn = 'work/W%d' %(scen)
        save(fn,W0)
        Hscen = GenHydroScen3d(nvec,nscen,nwells,scen)
    save('work/qhat',qhat)
    print('Setup Time=%4.4g secs') %(time()-t0)
#######################################################################
def phIteration3d(k,nscen,rho=-1,verbose=0):
   t0=time()
##   if rho==0:
##      cmd = 'mpiexec -n %d python scensolve3d.py 0' %(nscen)
##   else:
##      cmd = 'mpiexec -n %d python scensolve3d.py' %(nscen)
##   os.system(cmd)
   for scen in range(nscen):
       phScenSolve3d(scen)
   Ef=0;
   z = updateph3d()            # only uses nxyz if displaying Efk!!
   if verbose: print 'runtime[%d] = %4.4g secs' %(k,time()-t0)
   k+=1
   return k,Ef,z
#########################################################################
def knode3d(nx,ny,i,j,k):
    knode = nx*ny*k+ny*i+j
    return knode
##########################################################################
def HydroMatrix3d(nvec,nwells,baseKx,baseKy):
    nx,ny,nz=nvec;  nxyz=nx*ny*nz; #print 'NXYz=%d' %(nxyz)
    Kx = zeros((nxyz,1)); Ky = zeros((nxyz,1)); Kz = zeros((nxyz,1));
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                knode = knode3d(nx,ny,i,j,k)
                Kx[knode] = baseKx*(1.0+0.25/(1.0+i+j))
                Ky[knode] = baseKy*(1.0+0.25/(1.0+i+j))   # maybe replace by exp!
                Kz[knode] = (Kx[knode]+Ky[knode])/2.0
####################################################################
    data=[];ix=[];jx=[];
    # [1] Dirichlet in 2D: nodes 0->ny + ny*(nx-1)->nxy
    for i in range(nx*ny):
        k = i
        data.append(1);ix.append(k); jx.append(k); # Dirichlet UP!
##        data.append(1);ix.append(k); jx.append(k); # Dirichlet DN! should be Neu
    # (2) do all the rest of the nodes sumK*h-sum(K*h)=0
    for i in range(nx):
        for j in range(ny):
            for k in range(1,nz):
                knode = knode3d(nx,ny,i,j,k)
                Kdiag = 2*(Kx[knode]+Ky[knode]+Kz[knode])
                data.append(Kdiag);ix.append(knode); jx.append(knode); # Inner Diagonal
################################
                if (j>0):
                    data.append(-Kx[knode]);ix.append(knode); jx.append(knode-1); 
                else: #Neumann left
                    data.append(-Kx[knode]);ix.append(knode); jx.append(knode+1);                 
##############################
                if (i<nx-1):
                    data.append(-Ky[knode]);ix.append(knode); jx.append(knode+nx);
                else:
                    data.append(-Ky[knode]);ix.append(knode); jx.append(knode-nx);
####################
                if (k<nz-1):
                    data.append(-Kz[knode]);ix.append(knode); jx.append(knode+nx*ny);
                else:
                    data.append(-Kz[knode]);ix.append(knode); jx.append(knode-nx*ny);
####################
                if (j<ny-1):
                    data.append(-Kx[knode]);ix.append(knode); jx.append(knode+1);
                else: #Neumann right
                    data.append(-Kx[knode]);ix.append(knode); jx.append(knode-1);
                data.append(-Ky[knode]);ix.append(knode); jx.append(knode-nx);
####################
                data.append(-Kz[knode]);ix.append(knode); jx.append(knode-nx*ny);
####################                
    # (3) do the Neumann: fix by coo-additive prop: data[i][j]=-Kij
#    print 'done filling! NNZ=%d' %(len(data))
    A = coo_matrix((data,(ix,jx)),shape=(nxyz,nxyz))
    A = csc_matrix(A);
##    print A.todense()
##    print dot(A.todense(),ones((nxyz,1)))
    return A
####################################################################
def Presolve3d(nvec,nwells,baseKx,baseKy):
    t1=time();  nx,ny,nz=nvec; nxyz = nx*ny*nz;
    A = HydroMatrix3d(nvec,nwells,baseKx,baseKy)
    WellLocations = GenWellLocs3d(nvec,nwells)
############# set up RHS vector(s) ######################
    b=zeros((nxyz,1))
    for i in range(nx*ny):                   # check this well!
        b[i] = 0.4;
        b[nxyz-nx*ny+i] = 0.55
    b0 = zeros((nxyz,nwells))
############### solve! #####################
    xsol = zeros((nxyz,nwells+1))
    for i in range(nwells):
        b0[WellLocations[i],i]=1
        xsol[:,i] = linsolve.spsolve(A,b0[:,i],use_umfpack=False);
    xsol[:,nwells] = linsolve.spsolve(A,b,use_umfpack=False);
    return xsol
#########################################################################
def SolveOrLoadHwb3d(nx,ny,nz,scen,nwells,K0):
    fwb = 'Presolved/Hwb%d-%d-%d-%d-%d.npy' %(nx,ny,nz,nwells,scen)
    if os.path.exists(fwb):
       Hwb = load(fwb)
    else:           # else setup, solve, generate and save!
       baseKx=K0; baseKy=2*K0; baseKz = 0.5*K0
       Hwb = Presolve3d(nx,ny,nz,nwells,baseKx,baseKy,baseKz)
       save(fwb,Hwb)
    return Hwb
#########################################################################
def GenHydroScen3d(nvec,nscen,nwells,scen):
    WellLocations = GenWellLocs3d(nvec,nwells)
    Khyd = linspace(1.5,2.5,nscen)
    Hwb = SolveOrLoadHwb3d(nvec,scen,nwells,Khyd[scen])
    Dx = 0.011+0.01*scen;    Dy = 0.014+0.01*scen; Dz=0.0001+0.0001*scen;
    H = Hydrology3d(scen,Hwb,Dx,Dy,Dz,WellLocations)  
    if (scen==0):
        savetxt('work/WL',H.WL,'%d')
    return H
#########################################################################
def qhsave(q,h,scen):
    qname = 'work/q%d' %(scen);     savetxt(qname,q)
    hname = 'work/h%d' %(scen);     savetxt(hname,h)
#########################################################################
def headfrompump(q,H):
#    print 'in hfp',shape(H.Hwb),shape(q)
#    h = dot(H.Hwb[:,0:len(q)],q)+H.Hwb[:,-1]
#    print 'qhf:',q,shape(q)
    nwells = max(shape(q))
    h = H.Hwb[:,-1]
    qr= q.ravel()
#    print mean(h)
    for i in range(nwells):
        h = h + 10*qr[i]*H.Hwb[:,i]
#        print qr[i],mean(h)
#    print 'in hfp',shape(H.Hwb),shape(q),shape(h)
#    q=matrix(q);q=q.T
#    h = H.Hwb[:,0:len(q)]*q + H.Hwb[:,-1]
    return h
#######################################################################
def updateph3d(verbose=0):
    # updates qhat, W and maybe rho!
    rho,nscen,nwells,nx,ny,nz,LB,UB = load('work/refcard3d.npy');
#    print 'UPDATE: nscen =%d, nwells =%d' %(nscen,nwells)
    qhat = squeeze(zeros((nwells,1)));
    for i in range(nscen):
        qfname = 'work/qsol%d.npy' %(i);     
        qi = load(qfname);
##        print shape(qhat),shape(qi)
        qhat += qi/nscen
    save('work/qhat',qhat);
    if verbose: print 'qhat =',qhat.T
    znorm = 0
    for i in range(nscen):
        Wfname = 'work/W%d.npy' %(i);     Wi = squeeze(load(Wfname)); 
        qfname = 'work/qsol%d.npy' %(i);     qi = load(qfname); 
        Wi += rho*(qi-qhat);    save(Wfname,Wi)
        for j in range(nwells):
            znorm += (qi[j]-qhat[j])*(qi[j]-qhat[j])
    return znorm
############################################################################
def GenWellLocs3d(nvec,nwells):
    nx,ny,nz=nvec
    nx1q = 0.25*nx; nx2q = 0.5*nx; nx3q = 0.75*nx;  #assumes nx div by 4
    nzq = 0.75*nz # all wells at uniform depth!
    kpos = [(nx1q,nx1q),(nx1q,nx3q),(nx3q,nx1q),(nx3q,nx3q),
            (nx1q,nx2q),(nx2q,nx1q),(nx3q,nx2q),(nx2q,nx3q)]
    WellLocations = [knode3d(nx,ny,k[0],k[1],nzq) for k in kpos[:nwells]]   
    return WellLocations
#########################################################################
def nxyzload():
    rho,nscen,nwells,nx,ny,nz,LB,UB = load('work/refcard3d.npy');
    return nx,ny,nz    
#########################################################################
def GradientSum3d(h):
    sumgradx=0;sumgrady=0;sumgradz=0;
    nx,ny,nz = nxyzload()
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xij = knode3d(nx,ny,i,j,k);
                if i<nx-1:
                    gradhy_ij = h[xij+nx]-h[xij];
                    sumgrady = sumgrady + abs(gradhy_ij);
                if j<ny-1:
                    gradhx_ij = h[xij+1]-h[xij];
                    sumgradx = sumgradx + abs(gradhx_ij);
                if k<nz-1:
                    gradhz_ij = h[xij+nx*ny]-h[xij];
                    sumgradz = sumgradz + abs(gradhz_ij);

#    return 2*sumgradx + sumgrady
    return 10.0*(2*sumgradx + sumgrady + sumgradz)/(nx*ny*nz)
#######################################################################
def knode3d(nx,ny,i,j,k):
    return int(nx*ny*k+nx*i+j);
#########################################################################
def SolveOrLoadHwb3d(nvec,scen,nwells,K0):
    nx,ny,nz=nvec
    fwb = 'Presolved/Hwb%d-%d-%d-%d-%d.npy' %(nx,ny,nz,nwells,scen)
    if os.path.exists(fwb):
       Hwb = load(fwb)
    else:           # else setup, solve, generate and save!
       baseKx=K0; baseKy=2*K0;
       Hwb = Presolve3d(nvec,nwells,baseKx,baseKy)
       save(fwb,Hwb)
    return Hwb
#########################################################################
def phScenSolve3d(scen,verbose=0):
    Wfilename = 'work/W%d.npy' %(scen); W=load(Wfilename);
    qhat = load('work/qhat.npy');
    rho,nscen,nwells,nx,ny,nz,LB,UB = load('work/refcard3d.npy')
    H = GenHydroScen3d([nx,ny,nz],nscen,nwells,scen)
    e = lambda q,H,W,rho,qhat,scen: gwrfull(q,H,W,rho,qhat,scen)
    q0 = qhat; # we should read qsol!
    bounds = [(LB,UB) for i in range(len(qhat))]
    qopt = fmin_slsqp(e,q0,bounds=bounds,iprint=0,
                args=[H,W,rho,qhat,scen],acc=0.001)
    filename = 'work/qsol%d' %(scen);    save(filename,squeeze(qopt))
    print 'qsol[%d] =' %(scen),qopt.T
    if verbose: scenvecprint(scen,qopt)
    return
#########################################################################
def solplot3d(nvec,scen):
    nx,ny,nz=nvec
    use('Agg')
    import pylab
    fn = 'Presolved/Hwb%d-%d-%d-4-%d.npy' %(nx,ny,nz,scen)
    Hwb = load(fn)
    xi = Hwb[:,0]
    z=reshape(xi,nvec)
    for k in range(nz):
        zk = z[:,:,k]
        print zk
        pylab.subplot(3,3,k+1)
        pylab.imshow(z,interpolation='nearest')

    fn = 'FotoSlice-%d-%d-%d' %(nx,ny,nz)
    pyplot.savefig(fn)
########################################################################
def wrapitup3d(Ef,z,nvec,t1):
    k=len(Ef)
    print 'Converged in %d Iterations!' %(k)
    runtime = (time()-t1)/60
    print 'Ef(q%d)=%g    Solution Time=%4.2g minutes' %(k,Ef[-1],runtime)
#    ValidateSolution(25)
#    fn = 'Evo%d'  %(nx);savez(fn,z,Ef)
#    phplot(nx)
#########################################################################
def matshow(A):
    data = A.data
    nrows = len(A.indptr)
    for i in range(nrows):
        sta=A.indices[i]
        sto=A.indices[i+1]
        for j in range(sta,sto):
            print 'A[%d,%d]=%g' %(i,A.indptr[j],A.data[j]),
        print                                                    
#########################################################################
def gwrvec(q,H,scen):
#    print '*** invec',shape(q),shape(H.Hwb)
    h = headfrompump(q,H)
    qhsave(q,h,scen)  #necessary!
    hydraulic = GradientSum3d(h);
    perfcost = 0; #cpropagate2(q,h,H,scen);
#    pumping = sum(abs(q));
    pumping = sum(q);
    vec = [hydraulic, pumping, perfcost]
#    print 'hypype=',hydraulic,pumping,perfcost
    return vec
#########################################################################
def gwr(q,H,scen):
    " computes pumping costs and/or performance costs"
    vec = gwrvec(q,H,scen)
    y = sum(vec) #    y = hydraulic + pumping + perfcost
    return y
#########################################################################
def gwrfull(q,H,W,rho,qhat,scen):
    qdiff=(q-qhat).ravel()
    z = gwr(q,H,scen) + dot(q,W)+ 0.5*rho*dot(qdiff,qdiff)
#    print 'FIULLZ=',z
    return z
########################################
def GenHydrologies():
    H=[]
    rho,nscen,nwells,nx,ny,nz,LB,UB = load('work/refcard3d.npy');
    for i in range(nscen):
        H.append(GenHydroScen3d([nx,ny,nz],nscen,nwells,i))
    return H
#######################################################################
def hplot3d():
#    matplotlib.use('Agg')
#    import pylab
    qhat = load('work/qhat.npy')
    H = GenHydrologies()
    h = headfrompump(qhat,H[0])
    nx,ny,nz = nxyzload()
    z = reshape(h,(nx,ny,nz))
    save('z',z)
#    for layer in range(3):
#        pylab.subplot(3,3,layer+1)
#        pylab.imshow(z[layer],interpolation='nearest')
        
#    fn = 'FotoSurface';    matplotlib.pyplot.savefig(fn)

