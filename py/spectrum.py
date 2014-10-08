from numpy import *
from numpy.fft import fftn
npperp=16
npz=256
nproc=npperp*npz
i=0
nlx=1024
nly=1024
nlz=1024
nly_par=(nly-1)/npperp+1
nlz_par=(nlz-1)/npz+1
snapfold='/work/01022/tg802750/tod/1024/data_26_28'
snapname='snap0014.dat'

np=nlx/2
nkx=nlx
nky=nly
kx=zeros(nkx)
kx[:nkx/2+1]=arange(nkx/2+1)
kx[nkx/2+1:]=arange(-nkx/2+1,0)
ky=kx
zp=zeros((nlx,nly,nlz_par))
zm=zeros((nlx,nly,nlz_par))
sp=zeros(np)
sm=zeros(np)

for k in xrange(npz):
  print 'ipz='+str(k)
  for i in xrange(npperp):
    perpord=i%npperp
    dats=loadtxt(snapfold+'/'+str(i)+'/'+snapname)
    zp[:,perpord*nly_par:(perpord+1)*nly_par,:]=(dats[:,0].reshape(nlz_par,nly_par,nlx)).transpose()
    zm[:,perpord*nly_par:(perpord+1)*nly_par,:]=(dats[:,1].reshape(nlz_par,nly_par,nlx)).transpose()
  zpk=zeros((nkx,nky),dtype=complex)
  zmk=zeros((nkx,nky),dtype=complex)
  for zz in xrange(nlz_par):
    print 'spec ' +str(zz) + '/'+str(nlz_par)
    zpk=fftn(zp[:,:,zz])
    for ik in xrange(nkx):
      for jk in xrange(nky):
        kp=int(sqrt(kx[ik]**2+ky[jk]**2))
        if (kp<np):
          sp[kp]=sp[kp]+0.5*kp**2*(abs(zpk[ik,jk])**2)/(nkx*nky)**2/nlz
          sm[kp]=sm[kp]+0.5*kp**2*(abs(zmk[ik,jk])**2)/(nkx*nky)**2/nlz

spec=zeros((2,kp))    
spec[0,:]=sp
spec[1,:]=sm
savetxt(snapfold+'/spectrum_'+snapname+'.dat',spec)
