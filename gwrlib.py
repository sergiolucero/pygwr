from numpy import *
from scipy import sparse, linalg, io
from scipy.sparse import spdiags, csc_matrix, coo_matrix
from scipy.sparse.linalg import linsolve, bicg,cgs,bicgstab, gmres,lgmres
from scipy.optimize import fmin_l_bfgs_b, fmin_cobyla, fmin_slsqp
from openopt import NLP
from time import time
import os
import matplotlib
set_printoptions(suppress=True,precision=4)
########################################################################
class Hydrology(object):
    """ a hydrology class """
    def __init__ (self,n0scen,Hwb,Dx,Dy,WL):
        self.n0scen = n0scen; self.Hwb=Hwb;
        self.Dx=Dx; self.Dy=Dy; self.WL=WL
#########################################################################
class Refcard(object):
    """ a refcard class """  #too bad it can't be numpy-saved!
    def __init__ (self,rho,nscen,nwells,nx,LB,UB):
        self.rho = rho; self.nscen = nscen; self.nwells = nwells;
	self.nx = nx; self.LB = LB; self.UB = UB
#########################################################################
def knode(nx,i,j):
    return int(nx*i+j);
#########################################################################
def GradientSum(h):
    sumgradx=0;sumgrady=0;
    n2 = size(h); n=int(sqrt(n2));
    for i in range(n):
        for j in range(n):
            xij = knode(n,i,j);
            if i<n-1:
                gradhy_ij = h[xij+n]-h[xij];
                sumgrady = sumgrady + abs(gradhy_ij);
            if j<n-1:
                gradhx_ij = h[xij+1]-h[xij];
                sumgradx = sumgradx + abs(gradhx_ij);
    gsum = (2*sumgradx + sumgrady)/n2
#    print 'GS=%4.4g' %(gsum),
    return gsum
#########################################################################
def cpropagate(q,h,Hk):
#    Dx=Hk.Dx;Dy=Hk.Dy;
    nx = int(sqrt(size(h))); nwells = size(q)
#    try:
#        x= os.uname()		# should be done at the outset!!!
    cmd = './propagate %d %d %d' %(Hk.n0scen,nx,nwells)
#    except:
#        cmd = 'propagate %d %d %d' %(scen,nx,nwells)
    os.system(cmd)
    cfn = 'work/cost%d' %(Hk.n0scen)
    fp = open(cfn,"r")
    x = float(fp.read());
    fp.close()
#    print 'CP=',x
    return x   
#########################################################################
def sumsq(h):
    n=size(h)
    suma=0
    for i in range(n):
        suma += h[i]*h[i]
    return suma
#########################################################################
def headfrompump(q,H):
    nwells = size(q)
    h = H.Hwb[:,-1]
    qr= q.ravel()
    for i in range(nwells):
        h = h + 1000*qr[i]*H.Hwb[:,i]
#    print sum(h),sumsq(h)
    return h
#########################################################################
def gwrvec(q,H):
#    print '*** invec',shape(q),shape(H.Hwb)
    h = headfrompump(q,H)
    qhsave(q,h,H.n0scen)  #necessary for propagation!
    hydraulic = GradientSum(h);
    perfcost = cpropagate(q,h,H);
#    pumping = sum(abs(q));
    pumping = sum(q);
    vec = [pumping, hydraulic,perfcost]
    return vec
#########################################################################
def gwr(q,H):
    " computes pumping costs and/or performance costs"
#    wvec = [20,10.0,0.01*300]
    wvec = wvecload()
    gvec = gwrvec(q,H)
    y = dot(wvec,gvec)
#    y = sum(vec) #    y = hydraulic + pumping + perfcost
    return y
#########################################################################
def gwrCH(q,H):
    " computes pumping costs and/or performance costs"
#    wvec = [20,10.0,0.01*300]
    wvec = wvecload()
    gvec = gwrvec(q,H)
##    print 'inCHq[%d]=' %(scen),q
#    print 'q2=',sumsq(q)
#    print gvec
    y = dot(wvec,gvec)
#    y = sum(vec) #    y = hydraulic + pumping + perfcost
    return y
#########################################################################
def gwrfull(q,H,W,rho,qhat):
    qdiff=(q-qhat).ravel()
    z = gwr(q,H) + dot(q,W)+ 0.5*rho*dot(qdiff,qdiff)
    return z
############################################################################
def gwrfullCH(q,H,W,rho,qhat):
    qdiff=(q-qhat).ravel()
    z = gwr(q,H) + dot(q,W)+ 0.5*rho*dot(qdiff,qdiff)
    print 'GF=[%g,%g,%g]' %(gwrCH(q,H),dot(q,W),0.5*rho*dot(qdiff,qdiff))
    return z
############################################################################
def Expgwr(q,H):
    nscen = len(H)
    vec = zeros((nscen,3))
    for scen in range(nscen):
        vec[scen,:] = gwrvec(q,H[scen],scen)
#        suma += gwr(q,H[i])	# this should be done in parallel!
    suma = sum(vec[:])
    return suma/nscen
#########################################################################
def Expgwrverbose(q,H):
    nscen = len(H)
    vec = zeros((nscen,3))
    for scen in range(nscen):
        vec[scen,:] = gwrvec(q,H[scen],scen)
#        suma += gwr(q,H[i])	# this should be done in parallel!
    suma = sum(vec[:])
    print vec
    print suma/nscen
    return suma/nscen
#######################################################################
def GenerateWellLocations(nx,nwells):
    nx1q = 0.25*nx; nx2q = 0.5*nx; nx3q = 0.75*nx;  #assumes nx div by 4
    kpos = [(nx1q,nx1q),(nx1q,nx3q),(nx3q,nx1q),(nx3q,nx3q),
            (nx1q,nx2q),(nx2q,nx1q),(nx3q,nx2q),(nx2q,nx3q)]
    nwells = int(nwells); kpos = [(int(k[0]),int(k[1])) for k in kpos] #?!?!
    WellLocations = []
    
    for i in range(nwells):
	WellLocations.append(knode(nx,int(kpos[i][0]),int(kpos[i][1])))
#    WellLocations = [knode(nx,k[0],k[1]) for k in kpos[:nwells]]   
    return WellLocations
#######################################################################
def GenHydroScen(nx,nscen,nwells,scen):
    WellLocations = GenerateWellLocations(nx,nwells)
    Khyd = 1*linspace(2.5,3.5,nscen)
#    print 'K=',Khyd
    fwb = 'Presolved/Hwb%d-%d-%d.npy' %(nx,nwells,scen)
    if os.path.exists(fwb):
       Hwb = load(fwb)
    else:           # else setup, solve, generate and save!
       baseKx=Khyd[scen]; baseKy=2*Khyd[scen];
       Hwb = setupnsolve2d(nx,nx,nwells,baseKx,baseKy)
       save(fwb,Hwb)
    Dx = 0.011+0.01*scen;    Dy = 0.014+0.01*scen; 
    H = Hydrology(scen,Hwb,Dx,Dy,WellLocations)
    if (scen==0):
        savetxt('work/WL',H.WL,'%d')
    return H
#########################################################################
def GenHydrologies():
    H=[]
    rho,nscen,nwells,nx,LB,UB = readref()
    print type(nx)
    for i in range(nscen):
        H.append(GenHydroScen(nx,nscen,nwells,i))
    return H
#######################################################################
def scenvecprint(scen,q):
    q_ = q.reshape(1,-1)
    print 'qsol[%d] =[' %(scen),
    print q_.T
#########################################################################
def fmin_brute(q0,H,W,rho,qhat,scen):
    ntries = 100; radius = 0.01
    zmin = gwrfull(q0,H,W,rho,qhat)
    qbest = q0; qtest = q0
    nwells = max(shape(q0))
    print 'nw=%d' %(nwells)
    for i in range(ntries):
        qtest = qtest + radius*random.randn(nwells,1)
        ztest = gwrfull(qtest,H,W,rho,qhat)
        if (ztest<zmin):
            zmin = ztest
            qbest = qtest
    return qbest
#########################################################################
def readref():
    rho,nscen,nwells,nx,LB,UB = load('work/refcard.npy')
    nx = int(nx); nscen = int(nscen); nwells = int(nwells)
    return rho,nscen,nwells,nx,LB,UB
#########
def readcard(refcard):
    nx = refcard.nx
    LB = refcard.LB
    UB = refcard.UB
    return refcard.rho,refcard.nscen,refcard.nwells,nx,LB,UB
#########################################################################
def phScenSolve(scen,rho,verbose=0):
    Wfilename = 'work/W%d.npy' %(scen); W=load(Wfilename); 
    orho,nscen,nwells,nx,LB,UB = readref()
#    print 'rho=%g,nx=%d' %(rho,nx)
#    refcard = load('work/refcard.npy')
#    orho,nscen,nwells,nx,LB,UB = readcard(refcard)
    qfn = 'work/qhat%d.npy' %(nx); qhat = load(qfn)
    H = GenHydroScen(nx,nscen,nwells,scen)
    e = lambda q,H,W,rho,qhat: gwrfull(q,H,W,rho,qhat)
    q0 = qhat; # we should read qsol!
#    q0 = array([0.0025,-0.0038,0.0018,-0.0092])
    which_opt = 2 #0 slsqp,1 cobyla, 2 NLP
    if which_opt>0:
       if which_opt==1:
           up = lambda q,H,W,rho,qhat,scen: min(q-LB)
           dn = lambda q,H,W,rho,qhat,scen: max(UB-q)
           qopt = fmin_cobyla(e,q0,cons=[up,dn],iprint=0,
		    args=[H,W,rho,qhat],rhoend=0.0001)
#       qopt = fmin_brute(q0,H,W,rho,qhat,scen)
       else:
           eNLP = lambda q: gwrfull(q,H,W,rho,qhat)
           popt = NLP(eNLP,q0,lb=LB*ones((nwells,1)),
		ub=UB*ones((nwells,1)),iprint=-1)
           ropt = popt.solve('ralg')
           qopt = ropt.xf
           qopt = qopt.reshape(1,-1)
    else:
       bounds = [(LB,UB) for i in range(size(qhat))]
#       print bounds
       qopt = fmin_slsqp(e,q0,bounds=bounds,iprint=0,
		args=[H,W,rho,qhat,scen],acc=0.001)
   
    filename = 'work/qsol%d' %(scen);    save(filename,squeeze(qopt))
    print 'qsol[%d] =' %(scen),qopt
  #  qpert = zeros((1,nwells));
 #  for i in range(nwells):
##       qpert[:,i]= qopt[:,i]+0.01*random.randn()
#    print 'qpert[%d]=' %(scen),qpert
#    z1=gwrfullCH(qopt,H,W,rho,qhat,scen)
#    z2=gwrfullCH(qpert,H,W,rho,qhat,scen)
#    print 'TicToc=%g' %(z1-z2)
    if verbose: scenvecprint(scen,qopt)
    return
#######################################################################
def phIteration(k,nscen,rho=-1,verbose=0):
   t0=time()
   if rho==0:
      cmd = 'mpiexec -n %d python scensolve.py 0' %(nscen)
   else:
      cmd = 'mpiexec -n %d python scensolve.py' %(nscen)
   os.system(cmd)
   Ef,z = updateph()
   print 'z=',z
   if verbose: print 'runtime[%d] = %4.4g secs' %(k,time()-t0)
   k+=1
   return k,Ef,z
#######################################################################
def startph(nx,nscen,nwells,    LB = -1, UB = 1, rho = 10):
    t0=time()
    print '[Progressive Hedging Started][nx=%d,nscen=%d]' %(nx,nscen)
#    LB = -1; UB = 1; rho = 10;
    rho =1.0;
    W0=zeros((nwells,1));   qhat=matrix(zeros((nwells,1)))
#    refcard = Refcard(rho,nscen,nwells,nx,LB,UB)
#    save('work/refcard',refcard)
#    print refcard.nx
    save('work/refcard',[rho,nscen,nwells,nx,LB,UB])
    cmd = 'mpiexec -n %d python presolve.py' %(nscen)
    os.system(cmd)
    for scen in range(nscen):
        fn = 'work/W%d' %(scen);        save(fn,W0)
    qfn = 'work/qhat%d' %(nx);    save(qfn,qhat);
    print('Setup Time=%4.4g secs') %(time()-t0)
############################################################################
def updateph(verbose=0):
    # updates qhat, W and maybe rho!
    rho,nscen,nwells,nx,LB,UB = readref()
#    print 'UPDATE: nscen =%d, nwells =%d' %(nscen,nwells)
    qhat = zeros((1,nwells));
    for i in range(nscen):
        qfname = 'work/qsol%d.npy' %(i);     
        qi = load(qfname);    
        qhat += qi/nscen
    qfn = 'work/qhat%d' %(nx);    save(qfn,qhat);
    qhat=qhat.ravel()
    znorm = 0
    for i in range(nscen):
        Wfname = 'work/W%d.npy' %(i);     Wi = squeeze(load(Wfname)); 
        qfname = 'work/qsol%d.npy' %(i);     qi = load(qfname); 
        qi=qi.ravel();

        for j in range(nwells):
            Wi[j]+=rho*(qi[j]-qhat[j])
            znorm += (qi[j]-qhat[j])*(qi[j]-qhat[j])
            dq = qi[j]-qhat[j]
#	    print 'qij,qhj,diff,znorm',qi[j],qhat[j],dq,znorm
        save(Wfname,Wi)
    znorm = sqrt(znorm)
#   show progress (if wanted!)
#    H = GenHydrologies();
#    Ef = Expgwr(qhat,H)
#    print '>>> Efk = %4.6g zk = %2.3g' %(Ef,znorm)
    Ef=0
    return Ef,znorm
#########################################################################
def ValidateSolution(ntries = 10):
    t2=time()
    rho,nscen,nwells,nx,LB,UB = load('work/refcard.npy')
    qfn = 'work/qhat%d.npy' %(nx); qhat = load(qfn)
#    qhat=qhat.ravel(); qtest = qhat
    qhat=qhat.reshape(-1,1); qtest = qhat
    H = GenHydrologies();    
    fopt = Expgwrverbose(qhat,H);
    ftest = zeros((ntries,1))
   
    for i in range(ntries):
        for j in range(size(qhat)):
            qtest[j] = min(max(qhat[j] + 0.01*random.randn(),LB),UB)
        ftest[i] = Expgwrverbose(qtest,H)
        if ftest[i]<fopt:
            print 'OJO:',qtest.T
    fmin = min(ftest)
    nfailed = sum(ftest<fopt)
    print "%d failed out of %d, f(qopt) = %g, fmin=%g" %(nfailed,ntries,fopt,fmin)
    t3 = time();
    print "Validation Time =%2.2g mins" %((t3-t2)/60)
#########################################################################
def testhat():
    H=GenHydrologies()
    rho,nscen,nwells,nx,LB,UB = load('work/refcard.npy')
    qfn = 'work/qhat%d.npy' %(nx); qhat = load(qfn)
    Expgwrverbose(qhat,H)
#########################################################################
def wrapitup(Ef,z,nx,t1):
    k=len(Ef)
    print 'Converged in %d Iterations!' %(k)
    runtime = (time()-t1)/60
    print 'Ef(q%d)=%g    Solution Time=%4.2g minutes' %(k,Ef[-1],runtime)
    ValidateSolution(10)
#    fn = 'Evo%d'  %(nx);savez(fn,z,Ef)
#    phplot(nx)
#########################################################################
def qhsave(q,h,scen):
    qname = 'work/q%d' %(scen);     savetxt(qname,q)
    hname = 'work/h%d' %(scen);     savetxt(hname,h)
#########################################################################
def SystemMatrix(nx,ny,baseKx,baseKy):
    nxy = nx*ny;
    Kx = zeros((nxy,1));    Ky = zeros((nxy,1))
    for i in range(nx):
        for j in range(ny):
            k = ny*i + j
            Kx[k] = baseKx*(1.0+0.25/(1.0+i+j))
            Ky[k] = baseKy*(1.0+0.25/(1.0+i+j))   # maybe replace by exp!
    t0 = time()
    data=[];ix=[];jx=[];
    # [1] Dirichlet in 2D: nodes 0->ny + ny*(nx-1)->nxy
    for i in range(ny):
        k = i
        data.append(1);ix.append(k); jx.append(k); # Dirichlet UP!
        k = nxy-ny+i
        data.append(1);ix.append(k); jx.append(k); # Dirichlet DN!
    # (2) do all the rest of the nodes sumK*h-sum(K*h)=0
    for i in range(1,nx-1):
        for j in range(ny):
            k = knode(nx,i,j)
            Kdiag = 2*(Kx[k]+Ky[k])
            data.append(Kdiag);ix.append(k); jx.append(k); # Inner Diagonal
            if (j>0):
                data.append(-Kx[k]);ix.append(k); jx.append(k-1); 
            else: #Neumann left
                data.append(-Kx[k]);ix.append(k); jx.append(k+1);                 
            data.append(-Ky[k]);ix.append(k); jx.append(k+nx);
            if (j<ny-1):
                data.append(-Kx[k]);ix.append(k); jx.append(k+1);
            else: #Neumann right
                data.append(-Kx[k]);ix.append(k); jx.append(k-1);  
            data.append(-Ky[k]);ix.append(k); jx.append(k-nx);
            
    # (3) do the Neumann: fix by coo-additive prop: data[i][j]=-Kij

    A = coo_matrix((data,(ix,jx)),shape=(nxy,nxy))
    A = csc_matrix(A);

    return A
#########################################################################
def setupnsolve2d(nx,ny,nwells,baseKx,baseKy):
    nxy = nx*ny;
    A = SystemMatrix(nx,ny,baseKx,baseKy)
############# set up RHS vector(s) ######################
    WellLocations = GenerateWellLocations(nx,nwells)
    b=zeros((nxy,1))
    for i in range(nx):                   # check this well!
        b[i] = 0.4;        b[nx*(ny-1)+i] = 0.55
    b0 = zeros((nxy,nwells))
    b=zeros((nxy,1))
    for i in range(nx):
        b[i] = 0.4;        b[nx*(nx-1)+i] = 0.6
############### solve! #####################
    xsol = zeros((nxy,nwells+1))
    x0=ones((nxy,1))
    for i in range(nwells):
        b0[WellLocations[i],i]=1
        xsoli = LinSolve(A,b0[:,i],x0,'bicg')
#        print shape(xsol),shape(xsoli)
        xsol[:,i] = xsoli
    xsolnw = LinSolve(A,b,0.5*x0,'bicg');
    xsol[:,nwells] = xsolnw
##        xsol[:,i] = linsolve.spsolve(A,b0[:,i],use_umfpack=True);
##    xsol[:,nwells] = linsolve.spsolve(A,b,use_umfpack=True);
##    xplot(xsol)
    return xsol
######## PLOTTING ROUTINES ############################################
def LinSolve(A,b,x0,method):
    if method=='bicg':
        x,info = bicg(A,b,x0)
    if method=='cgs':
        x,info = cgs(A,b,x0)
    if method=='bicgstab':
        x,info = bicgstab(A,b,x0)
        if (info): print 'INFO=',info
    if method=='superLU':
        x = linsolve.spsolve(A,b,use_umfpack=False) 
    if method=='umfpack':
        x = linsolve.spsolve(A,b,use_umfpack=True) 
    if method=='gmres':
        x,info = gmres(A,b,x0) 
    if method=='lgmres':
        x,info = lgmres(A,b,x0) 
    return x
########################################
def solplot(nx,nwells,scen):
    matplotlib.use('Agg')
    import pylab
    fn = 'Presolved/Hwb%d-4-%d.npy' %(nx,scen)
    Hwb = load(fn)
    for i in range(nwells+1):
        xi = Hwb[:,i]
        pylab.subplot(3,3,i+1)
        z=reshape(xi,(nx,nx))
        pylab.imshow(z,interpolation='nearest')
        fn = 'FotoWell%d-%d-%d' %(nx,nwells,scen)
        matplotlib.pyplot.savefig(fn)
############################################################
def phplot(nx):
    matplotlib.use('Agg');
    import pylab
    fn = 'Evo%d.npz' %(nx)
    x=load(fn);Ef=x['arr_1'];z=x['arr_0']
    pylab.subplot(211);pylab.plot(z)
    pylab.subplot(212);pylab.plot(Ef)
    fn2 = 'Foto%d.png' %(nx);
    pylab.savefig(fn2)
#########################################################################
def hplot():
    matplotlib.use('Agg');    import pylab
    x=load('hlast.npy');
    nx2,ny=shape(x)
    nx=int(sqrt(nx2))
    x0=reshape(x[0,:],(nx,nx))
    pylab.plot(x0)
    fn2 = 'hFoto%d.png' %(nx);
    pylab.savefig(fn2)
#########################################################################
def nxload():
    rho,nscen,nwells,nx,LB,UB = load('work/refcard.npy')
    return nx
#########################################################################
def wvecload():
    w=load('work/wvec.npy')
    return w
#########################################################################
def qhatload():
    nx = nxload()
    qfn = 'work/qhat%d.npy' %(nx); 
    qhat = load(qfn)
    return qhat
#########################################################################
def SolTest():
    qhat = qhatload()
    H=GenHydrologies()
    Expgwrverbose(qhat,H)
    for i in range(5):
        qhat = qhat + random.randn(size(qhat))
        print qhat
        Expgwrverbose(qhat,H)
    return
#########################################################################
def OneRun(nx):
    nscen = 4
    cmd = 'python parph.py %d %d 4' %(nx,nscen)
    os.system(cmd)
#########################################################################
def FullRun():
    for i in range(8,12):
        nx = 2**i
        OneRun(nx)
#########################################################################
def LinTest(baseKx=1.0):
    nwells = 4;
    baseKy = 2.0*baseKx
    narray=[2**i for i in range(5,9)]
    methods = ['bicg','bicgstab','cgs','gmres','lgmres','superLU','umfpack']
#    narray=[100*i for i in range(3,11)]
#    methods = ['bicg','bicgstab','superLU','umfpack']
    nd=len(narray);nm=len(methods);rez=zeros((nd,nm))

    for j in range(nd):
        nx=narray[j]
        print nx
        A = SystemMatrix(nx,nx,baseKx,baseKy)
        b = random.randn(nx*nx,nwells+1)
#        x0 = zeros((nx*nx,1))   # COMPARE AGAINST X0 = AVG(B)*SIZE(B)
        x0 = mean(b)*ones((nx*nx,1))   # COMPARE AGAINST X0 = AVG(B)*SIZE(B)
        for i in range(nm):
            t1 = time()
            for k in range(nwells+1):
                x = LinSolve(A,b[:,k],x0,methods[i])
#            print '[%s] %4.2g ' %(method, time()-t1),
            rez[j,i]=(time()-t1)
    return rez
########################################################################        
def FullTest():
    ntries = 5
    rez = LinTest()
    for k in range(ntries-1):
        rez += LinTest()
    return rez/ntries
########################################################################        
def mp(A):
    m,n = shape(A)
    for i in range(m):
        print 2**(i+5),' & ',
        for j in range(n):
            print '%4.2g' %(A[i,j]),
            print ' & ',
        print '\\'
########################################################################        
def make():
    cmd = 'make'
    os.system(cmd)
########################################################################        
def run(n):
    cmd = 'python parph.py %d 4 4' %(n)
    os.system(cmd)
########################################################################        
def cedit():
    cmd = 'nano prop.c'
    os.system(cmd)
########################################################################        
#def relo():
#    from gwrlib import *
#################################
def qh():
    q=qhatload()
    H = GenHydrologies()
    return q,H
################
def convextest(npoints=10):
    q,H = qh()
    qmin = -1.0*ones((4,1))
    qmax = 0.7*ones((4,1))
    alfa = linspace(0,1,npoints)

    for i in range(npoints):
        qtest = (1-alfa[i])*qmin + alfa[i]*qmax
        ztest = Expgwr(qtest,H)
        print '%2.2g,%4.4g' %(alfa[i],ztest)
###########################################################################
def noway():
    qhat,H = qh()
    ntries = 100
    nwells = max(shape(qhat))
    print 'QHAT',qhat,Expgwr(qhat,H)

    for i in range(ntries):
        qtest = (-1.0 + 2.0*random.rand(nwells,1))/10.0
        z = Expgwr(qtest,H)
        if (z<0.1):
            print '%2.2g,%4.2g' %(qtest.T, z)
###########################################################################
def cmovie(q,scen):
    rho,nscen,nwells,nx,LB,UB = load('work/refcard.npy')
    cmd = './plotpropagate %d %d %d' %(scen,nx,nwells)
    print cmd
    os.system(cmd)
###########################################################################
def softsetup():
   os.system('make -s')  # updates propagator!
   os.system('easy_install --quiet mpi4py'); # if not already available!
   os.system('apt-get -qq install zip');
   os.system('svn -q co svn://openopt.org/PythonPackages OOSuite')
   opendir = '/usr/local/lib/python2.6/dist-packages/openopt-0.31-py2.6.egg/'
   if not(os.path.isdir(opendir)): os.system('cd OOSuite;python install_all.py')
   if not(os.path.isdir('work')): os.mkdir('work');
   if not(os.path.isdir('Presolved')): os.mkdir('Presolved');

