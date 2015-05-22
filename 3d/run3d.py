import gwrlib3d

nx=32;ny=32;nz=10;
nscen=4; nwells=4
print 'full-dims = %d' %(nx*ny*nz)

gwrlib3d.startph3d([nx,ny,nz],nscen,nwells)		# now done in parallel
k,Efk,zk = gwrlib3d.phIteration3d(0,nscen,0,0)  # k=0, rho=0
#gwrlib3d.wrapitup3d
