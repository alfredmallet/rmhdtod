from os import listdir
mode=0
mode2=0
nm=10
imin=10
iminperp=10
imaxperp=30
imax=30
iminpar=10
imaxpar=27
ddir='1024_sfdata'
#ddir='1024_sfdata_neg'
lsd=listdir(ddir)
nsn=len(lsd)
isn=0
seps=zeros(32)
avsf=zeros((32,9,4,9,nm))
avsf2=zeros(shape(avsf))
nsftot=zeros((32,9,4,9,nm))
astoreperp=zeros((nsn,nm))
astorefluc=zeros((nsn,nm))
astorepar=zeros((nsn,nm))
astorefull=zeros((nsn,nm))
avpar2=zeros((32,nm))
avpar3=zeros((32,nm))
for isep in xrange(32):
  seps[isep]=4.*(1024.*0.25/4)**((isep)*1./31.)
if mode==1:
    imin=7
    iminperp=7
    imax=27
    imaxperp=27
for lss in lsd:
  print isn, lss
  d=loadtxt(ddir+'/'+lss)
  zpsf=d[:,6+mode]
  print zpsf[0]
  nsf=d[:,5]
  nsfr=nsf.reshape((32,9,4,9,nm))
  zpsfr=zpsf.reshape((32,9,4,9,nm))
  nzpsfr=(nsf*zpsf).reshape((32,9,4,9,nm))
  avsf=avsf+nzpsfr
  avsf2=avsf2+(1./nsn)*log(zpsfr)
  nsftot=nsftot+nsfr
  lam=seps*pi/512.
  perp=zpsfr[:,-1,mode2,-1,:]
  fluc=zpsfr[:,-1,mode2,0,:]
  par=zpsfr[:,0,mode2,-1,:]
  full=sum(nzpsfr[:,-1,mode2,:,:],axis=1)/sum(nsfr[:,-1,mode2,:,:]*1.0,axis=1)
  figure(1)
  loglog(lam,par[:,-1])
  avpar2=avpar2+par/nsn
  avpar3=avpar3+(1./nsn)*log(par)
  #full=sum(sum(nzpsfr[:,:,mode,:,:],axis=1),axis=1)/sum(sum(nsfr[:,:,mode,:,:]*1.0,axis=1),axis=1)
  for m in xrange(0,nm):
    print m
    a=polyfit(log(lam[iminperp:imaxperp]),log(perp[iminperp:imaxperp,m]),1)[0]
    astoreperp[isn,m]=a
    a=polyfit(log(seps[imin:imax]),log(fluc[imin:imax,m]),1)[0]
    astorefluc[isn,m]=a
    a=polyfit(log(seps[iminpar:imaxpar]),log(par[iminpar:imaxpar,m]),1)[0]
    astorepar[isn,m]=a
    a=polyfit(log(seps[iminperp:imaxperp]),log(full[iminperp:imaxperp,m]),1)[0]
    astorefull[isn,m]=a
  isn=isn+1

avpar3=exp(avpar3)
avsf=avsf/nsftot
avsf2=exp(avsf2)
avfull=sum(avsf[:,-1,mode,:,:]*nsftot[:,-1,mode,:,:],axis=1)/sum(1.0*nsftot[:,-1,mode,:,:],axis=1)
avperp=avsf2[:,-1,mode,-1,:]
avfluc=avsf2[:,-1,mode,0,:]
avpar=avsf2[:,0,mode,-1,:]
aperps=zeros(nm)
aflucs=zeros(nm)
apars=zeros(nm)
afulls=zeros(nm)
afull=zeros(nm)
siperps=zeros(nm)
siflucs=zeros(nm)
sipars=zeros(nm)
sifulls=zeros(nm)
for m in xrange(0,nm):
  aperps[m]=polyfit(log(lam[iminperp:imaxperp]),log(avperp[iminperp:imaxperp,m]),1)[0]
  aflucs[m]=polyfit(log(seps[imin:imax]),log(avfluc[imin:imax,m]),1)[0]
  apars[m]=polyfit(log(seps[iminpar:imaxpar]),log(avpar[iminpar:imaxpar,m]),1)[0]
  afulls[m]=polyfit(log(seps[imin:imax]),log(avfull[imin:imax,m]),1)[0]

#aperps=mean(astoreperp,axis=0)
#aflucs=mean(astorefluc,axis=0)
#apars=mean(astorepar,axis=0)
#afulls=mean(astorefull,axis=0)
errfact=2./sqrt(nsn)
#errfact=1.
siperps2=zeros(nm)
siflucs2=zeros(nm)
sipars2=zeros(nm)
ddlperp=((log(avperp[:-1,:])-log(avperp[1:,:])).transpose()/(log(seps[:-1])-log(seps[1:]))).transpose()
ddlfluc=((log(avfluc[:-1,:])-log(avfluc[1:,:])).transpose()/(log(seps[:-1])-log(seps[1:]))).transpose()
ddlpar=((log(avpar[:-1,:])-log(avpar[1:,:])).transpose()/(log(seps[:-1])-log(seps[1:]))).transpose()
siperps3=zeros(nm)
siflucs3=zeros(nm)
sipars3=zeros(nm)
for m in xrange(nm):
  siperps[m]=errfact*sqrt(sum((aperps[m]-astoreperp[:,m])**2)/(1.0*nsn))
  siflucs[m]=errfact*sqrt(sum((aflucs[m]-astorefluc[:,m])**2)/(1.0*nsn))
  sipars[m]=errfact*sqrt(sum((apars[m]-astorepar[:,m])**2)/(1.0*nsn))
  sifulls[m]=errfact*sqrt(sum((afulls[m]-astorefull[:,m])**2)/(1.0*nsn))
  siperps2[m]=(polyfit(log(lam[iminperp:(iminperp+imaxperp)/2]),log(avperp[iminperp:(iminperp+imaxperp)/2,m]),1)[0]-polyfit(log(lam[(iminperp+imaxperp)/2:imaxperp]),log(avperp[(iminperp+imaxperp)/2:imaxperp,m]),1)[0])/2
  siflucs2[m]=(polyfit(log(lam[imin:(imin+imax)/2]),log(avfluc[imin:(imin+imax)/2,m]),1)[0]-polyfit(log(lam[(imin+imax)/2:imax]),log(avfluc[(imin+imax)/2:imax,m]),1)[0])/2
  sipars2[m]=(polyfit(log(lam[iminpar:(iminpar+imaxpar)/2]),log(avpar[iminpar:(iminpar+imaxpar)/2,m]),1)[0]-polyfit(log(lam[(iminpar+imaxpar)/2:imaxpar]),log(avpar[(iminpar+imaxpar)/2:imaxpar,m]),1)[0])/2
  siperps3[m]=sqrt(mean((ddlperp[iminperp:imaxperp,m]-aperps[m])**2))/sqrt(imaxperp-iminperp*1.0)
  siflucs3[m]=sqrt(mean((ddlfluc[imin:imax,m]-aflucs[m])**2))/sqrt(imax-imin*1.0)
  sipars3[m]=sqrt(mean((ddlpar[iminpar:imaxpar,m]-apars[m])**2))/sqrt(imaxpar-iminpar*1.0)
figure()
loglog(lam,avperp)
for m in xrange(nm):
  loglog(lam,(avperp[imax,m]/lam[imax]**aperps[m])*lam**aperps[m],'k:')
figure()
loglog(lam,avfluc)
for m in xrange(nm):
  loglog(lam,(avfluc[imax,m]/lam[imax]**aflucs[m])*lam**aflucs[m],'k:')
figure()
loglog(lam,avpar)
for m in xrange(nm):
  loglog(lam,(avpar[imax,m]/lam[imax]**apars[m])*lam**apars[m],'k:')
figure()
loglog(lam,avfull)
for m in xrange(nm):
    loglog(lam,(avfull[imax,m]/lam[imax]**afulls[m])*lam**afulls[m],'k:')
mstore=arange(0.0,nm/2.+0.5,.5)
ms=mstore[1:]
mstore2=arange(0.0,10.0,.5)
figure(figsize=array([9.5,8]))
#plot(mstore2, 1-0.691**mstore2,'k:')
#plot(mstore2,1-sqrt(0.5)**mstore2,'k:')
#plot(mstore2,mstore2*(1-sqrt(0.5)**mstore2)/(mstore2/2+1-sqrt(0.5)**mstore2),'k:')
#plot(mstore2,2*(1-sqrt(0.5)**mstore2),'k:')
errorbar(mstore,concatenate([[0],aperps]),yerr=concatenate([[0],siperps]),color='b',linewidth=1.5)
errorbar(mstore,concatenate([[0],aflucs]),yerr=concatenate([[0],siflucs]),color='g',linewidth=1.5)
errorbar(mstore,concatenate([[0],apars]),yerr=concatenate([[0],sipars]),color='r',linewidth=1.5)
#errorbar(mstore,concatenate([[0],afulls]),yerr=concatenate([[0],sifulls]),color='k',linewidth=1.5)
beta=0.691
kf=mean(aflucs/(1-beta**ms))
kpar=mean(apars/(1-beta**ms))
#kf=1.41
#kpar=1/(1-beta**2)
print 'kf ', kf
print 'kpar', kpar
plot(mstore2,1-beta**mstore2,'k:')
plot(mstore2,kpar*(1-beta**mstore2),'k:')
#plot(mstore2,2*(1-beta**mstore2)/(1-2*beta**mstore2 * log(beta)),'k:')
plot(mstore2,kf*(1-beta**mstore2),'k:')
xlim(0,nm/2.+0.5)
#text(0.5,1.75,r'Total',color='k',fontsize=20)
text(3,0.5,r'$\zeta^\perp_n$',color='b',fontsize=20)
text(4.07,1,r'$\zeta^{fluc}_n$',color='g',fontsize=20)
text(2.5,1.4,r'$\zeta^\parallel_n$',color='r',fontsize=20)
xlabel(r'$n$',fontsize=30)
ylabel(r'$\zeta_n$',fontsize=30)
#figure()
show()


