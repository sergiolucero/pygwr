import matplotlib
from numpy import *
matplotlib.use('Agg')
from gwrlib import solplot, seepropagate, GenHydrologies
from gwrlib import qhsave, headfrompump, cpropagate2, nxload
from matplotlib.backends.backend_pdf import PdfPages
##################################################
def m1():
    pp = PdfPages('multipage.pdf')

    nx=32
    for i in range(3):
        solplot(nx,4,0);
        titulo = 'Nx=%d' %(nx);   matplotlib.pyplot.title(titulo);
        pp.savefig()
        nx=nx*2
    pp.close()
##################################################
def m2():
    pp = PdfPages('m2.pdf')
    import pylab
    nsteps=10
    qhat = load('work/qhat.npy'); qtest = qhat
    H = GenHydrologies();    h = headfrompump(qhat,H[0])
    nx = int(sqrt(len(h)));  nc = nx*(nx/2-1)+nx/2;
    spot = nc+array([-nx,-nx+1,-1,0,1,nx-1,nx]);
    cprof = zeros((nx*nx,1));    cnew = zeros((nx*nx,1));
    for j in spot: cprof[j]=0.5;
    for i in range(nsteps):
        cprof = seepropagate(cprof,i,qhat,h,H[0])
        fn = 'fotos/Cprof%d.npy' %(i); ci = load(fn)
        z = reshape(ci,(nx,nx))
        pylab.imshow(z,interpolation='nearest')
        pp.savefig()
#        matplotlib.pyplot.title(titulo);
    pp.close()
##################################################
def m3():
    qhat = load('work/qhat32.npy'); 
    dum,nwells = shape(qhat)
    H = GenHydrologies();    nx2,nwp = shape(H[0].Hwb)
    nx = int(sqrt(nx2))

    ntries = 10
    for i in range(ntries):
        qtest = array(random.randn(nwells,1))
        qtest = -0.5+0.1*i*ones((nwells,1))
        qtest=qtest.T
        h = headfrompump(qtest,H[0])
	qhsave(qtest,h,0)
        ctest = cpropagate2(qtest,h,H[0],0)
        print qtest,'cT(q)=',ctest
    save('hlast',h)
#####################
m3()
