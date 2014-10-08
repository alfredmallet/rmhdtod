from os import listdir

mode=0
imin=10
iminperp=10
imax=27
iminpar=10
imaxpar=27
lsd=listdir('512_sfdata')
nsn=len(lsd)
isn=0
seps=zeros(32)
avsf=zeros((32,9,4,9,10))
nsftot=zeros((32,9,4,9,10))
astoreperp=zeros((nsn,10))
astorefluc=zeros((nsn,10))
astorepar=zeros((nsn,10))
astorefull=zeros((nsn,10))
for isep in xrange(32):
  seps[isep]=4.*(512.*0.25/4)**((isep)*1./31.)
for lss in lsd:
  print isn, lss
  d=loadtxt('512_sfdata/'+lss+'/3dsf.dat')
  zpsf=d[:,6+mode]
  nsf=d[:,5]
  nsfr=nsf.reshape((32,9,4,9,10))
  zpsfr=zpsf.reshape((32,9,4,9,10))
  nzpsfr=(nsf*zpsf).reshape((32,9,4,9,10))
  avsf=avsf+nzpsfr
  nsftot=nsftot+nsfr
  lam=seps*pi/256.
  perp=zpsfr[:,-1,mode,-1,:]
  fluc=zpsfr[:,-1,mode,0,:]
  par=zpsfr[:,0,mode,-1,:]
  #full=sum(nzpsfr[:,-1,mode,:,:],axis=1)/sum(nsfr[:,-1,mode,:,:]*1.0,axis=1)
  full=sum(sum(nzpsfr[:,:,mode,:,:],axis=1),axis=1)/sum(sum(nsfr[:,:,mode,:,:]*1.0,axis=1),axis=1)
  for m in xrange(0,10):
    print m
    a=polyfit(log(lam[iminperp:imax]),log(perp[iminperp:imax,m]),1)[0]
    astoreperp[isn,m]=a
    a=polyfit(log(seps[imin:imax]),log(fluc[imin:imax,m]),1)[0]
    astorefluc[isn,m]=a
    a=polyfit(log(seps[iminpar:imaxpar]),log(par[iminpar:imaxpar,m]),1)[0]
    astorepar[isn,m]=a
    a=polyfit(log(seps[imin:imax]),log(full[imin:imax,m]),1)[0]
    astorefull[isn,m]=a
  isn=isn+1

avsf=avsf/nsn
avfull=sum(avsf[:,-1,mode,:,:],axis=1)/sum(1.0*nsftot[:,-1,mode,:,:],axis=1)
avperp=avsf[:,-1,mode,-1,:]
avfluc=avsf[:,-1,mode,0,:]
avpar=avsf[:,0,mode,-1,:]
aperps=zeros(10)
aflucs=zeros(10)
apars=zeros(10)
afulls=zeros(10)
afull=zeros(10)
siperps=zeros(10)
siflucs=zeros(10)
sipars=zeros(10)
sifulls=zeros(10)
for m in xrange(0,10):
  aperps[m]=polyfit(log(lam[imin:imax]),log(avperp[imin:imax,m]),1)[0]
  aflucs[m]=polyfit(log(seps[imin:imax]),log(avfluc[imin:imax,m]),1)[0]
  apars[m]=polyfit(log(seps[iminpar:imaxpar]),log(avpar[iminpar:imaxpar,m]),1)[0]
  afulls[m]=polyfit(log(seps[imin:imax]),log(avfull[imin:imax,m]),1)[0]

aperps=mean(astoreperp,axis=0)
aflucs=mean(astorefluc,axis=0)
apars=mean(astorepar,axis=0)
afulls=mean(astorefull,axis=0)
for m in xrange(10):
  siperps[m]=sqrt(sum((aperps[m]-astoreperp[:,m])**2)/(1.0*nsn))
  siflucs[m]=sqrt(sum((aflucs[m]-astorefluc[:,m])**2)/(1.0*nsn))
  sipars[m]=sqrt(sum((apars[m]-astorepar[:,m])**2)/(1.0*nsn))
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
loglog(lam,avpar)
for m in xrange(10):
  loglog(lam,(avpar[imax,m]/lam[imax]**apars[m])*lam**apars[m],'k:')
figure()
loglog(lam,avfull)
for m in xrange(10):
    loglog(lam,(avfull[imax,m]/lam[imax]**afulls[m])*lam**afulls[m],'k:')
mstore=arange(0.0,5.5,.5)
mstore2=arange(0.0,10.0,.5)
figure()
#plot(mstore2, 1-0.691**mstore2,'k:')
plot(mstore2,1-sqrt(0.5)**mstore2,'k:')
plot(mstore2,mstore2*(1-sqrt(0.5)**mstore2)/(mstore2/2+1-sqrt(0.5)**mstore2),'k:')
plot(mstore2,2*(1-sqrt(0.5)**mstore2),'k:')
errorbar(mstore,concatenate([[0],aperps]),yerr=concatenate([[0],siperps]),color='b',linewidth=1.5)
errorbar(mstore,concatenate([[0],aflucs]),yerr=concatenate([[0],siflucs]),color='g',linewidth=1.5)
errorbar(mstore,concatenate([[0],apars]),yerr=concatenate([[0],sipars]),color='r',linewidth=1.5)
errorbar(mstore,concatenate([[0],afulls]),yerr=concatenate([[0],sifulls]),color='k',linewidth=1.5)
xlim(0,5.5)
text(0.5,1.75,r'Total',color='k',fontsize=20)
text(0.5,1.6,r'$\perp$',color='b',fontsize=20)
text(0.5,1.45,r'Fluctuation',color='g',fontsize=20)
text(0.5,1.3,r'$\parallel$',color='r',fontsize=20)
show()
