from numpy import *
from pylab import *
from os import listdir

mode=0
imin=20
iminperp=15
imax=31
iminpar=15
imaxpar=30
lsd=listdir('2d_512_sfdata')
nsn=len(lsd)
isn=0
seps=zeros(32)
avsf=zeros((32,4,9,10))
nsftot=zeros((32,4,9,10))
astoreperp=zeros((nsn,10))
astorefluc=zeros((nsn,10))
astorepar=zeros((nsn,10))
astorefull=zeros((nsn,10))
for isep in xrange(32):
  seps[isep]=4.*(512.*0.25/4)**((isep)*1./31.)
for lss in lsd:
  print isn, lss
  d=loadtxt('2d_512_sfdata/'+lss)
  zpsf=d[:,5+mode]
  nsf=d[:,4]
  nsfr=nsf.reshape((32,4,9,10))
  zpsfr=zpsf.reshape((32,4,9,10))
  nzpsfr=(nsf*zpsf).reshape((32,4,9,10))
  avsf=avsf+nzpsfr
  nsftot=nsftot+nsfr
  lam=seps*pi/256.
  perp=zpsfr[:,mode,-1,:]
  fluc=zpsfr[:,mode,0,:]
  full=(sum(nzpsfr[:,mode,:,:],axis=1))/(sum(nsfr[:,mode,:,:]*1.0,axis=1))
  for m in xrange(0,10):
    print m
    a=polyfit(log(lam[iminperp:imax]),log(perp[iminperp:imax,m]),1)[0]
    astoreperp[isn,m]=a
    a=polyfit(log(seps[imin:imax]),log(fluc[imin:imax,m]),1)[0]
    astorefluc[isn,m]=a
    a=polyfit(log(seps[imin:imax]),log(full[imin:imax,m]),1)[0]
    astorefull[isn,m]=a
  isn=isn+1
  #figure(100)
  #loglog(lam,full[:,3])

avsf=avsf/nsftot
avfull=sum(avsf[:,mode,:,:]*nsftot[:,mode,:,:],axis=1)/sum(1.0*nsftot[:,mode,:,:],axis=1)
avperp=avsf[:,mode,-1,:]
avfluc=avsf[:,mode,0,:]
aperps=zeros(10)
aflucs=zeros(10)
afulls=zeros(10)
afull=zeros(10)
siperps=zeros(10)
siflucs=zeros(10)
sifulls=zeros(10)
for m in xrange(0,10):
  aperps[m]=polyfit(log(lam[imin:imax]),log(avperp[imin:imax,m]),1)[0]
  aflucs[m]=polyfit(log(seps[imin:imax]),log(avfluc[imin:imax,m]),1)[0]
  afulls[m]=polyfit(log(seps[imin:imax]),log(avfull[imin:imax,m]),1)[0]

aperps=mean(astoreperp,axis=0)
aflucs=mean(astorefluc,axis=0)
afulls=mean(astorefull,axis=0)
for m in xrange(10):
  siperps[m]=sqrt(sum((aperps[m]-astoreperp[:,m])**2)/(1.0*nsn))
  siflucs[m]=sqrt(sum((aflucs[m]-astorefluc[:,m])**2)/(1.0*nsn))
  sifulls[m]=sqrt(sum((afulls[m]-astorefull[:,m])**2)/(1.0*nsn))


figure()
loglog(lam,avperp)
for m in xrange(10):
  loglog(lam,(avperp[imax,m]/lam[imax]**aperps[m])*lam**aperps[m],'k:')
figure()
loglog(lam,avfluc)
for m in xrange(10):
  loglog(lam,(avfluc[imax,m]/lam[imax]**aflucs[m])*lam**aflucs[m],'k:')
figure()
loglog(lam,avfull)
for m in xrange(10):
    loglog(lam,(avfull[imax,m]/lam[imax]**afulls[m])*lam**afulls[m],'k:')
mstore=arange(0.0,5.5,.5)
mstore2=arange(0.0,10.0,.5)
figure()
plot(mstore2, 1-0.691**mstore2,'k:')
#plot(mstore2,1-sqrt(0.5)**mstore2,'k:')
plot(mstore2,mstore2*(1-sqrt(0.5)**mstore2)/(mstore2/2+1-sqrt(0.5)**mstore2),'k:')
errorbar(mstore,concatenate([[0],aperps]),yerr=concatenate([[0],siperps]),color='b',linewidth=1.5)
errorbar(mstore,concatenate([[0],aflucs]),yerr=concatenate([[0],siflucs]),color='g',linewidth=1.5)
errorbar(mstore,concatenate([[0],afulls]),yerr=concatenate([[0],sifulls]),color='k',linewidth=1.5)
xlim(0,5.5)
text(0.5,1.75,r'Total',color='k',fontsize=20)
text(0.5,1.6,r'$\perp$',color='b',fontsize=20)
text(0.5,1.45,r'Fluctuation',color='g',fontsize=20)
show()
