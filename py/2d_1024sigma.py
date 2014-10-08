from os import listdir

mode=0
imin=10
print imin
iminperp=imin
imax=30
sfdir='2dsfdata_m20'
#sfdir='2d_1024_sfdata'
lsd=listdir(sfdir)
nsn=len(lsd)
isn=0
nm=20
seps=zeros(32)
avsf=zeros((32,4,9,nm))
nsftot=zeros((32,4,9,nm))
astoreperp=zeros((nsn,nm))
astorefluc=zeros((nsn,nm))
astorepar=zeros((nsn,nm))
astorefull=zeros((nsn,nm))
for isep in xrange(32):
  seps[isep]=4.*(1024.*0.25/4)**((isep)*1./31.)
for lss in lsd[:-1]:
  #print isn, lss
  d=loadtxt(sfdir+'/'+lss)
  zpsf=d[:,5+mode]
  nsf=d[:,4]
  nsfr=nsf.reshape((32,4,9,nm))
  zpsfr=zpsf.reshape((32,4,9,nm))
  nzpsfr=(nsf*zpsf).reshape((32,4,9,nm))
  #avsf=avsf+nzpsfr
  avsf=avsf+log(zpsfr)/nsn
  nsftot=nsftot+nsfr
  lam=seps*pi/512.
  perp=zpsfr[:,mode,-1,:]
  fluc=zpsfr[:,mode,0,:]
  full=(sum(nzpsfr[:,mode,:,:],axis=1))/(sum(nsfr[:,mode,:,:]*1.0,axis=1))
  for m in xrange(0,nm):
    #print m
    a=polyfit(log(lam[iminperp:imax]),log(perp[iminperp:imax,m]),1)[0]
    astoreperp[isn,m]=a
    a=polyfit(log(seps[imin:imax]),log(fluc[imin:imax,m]),1)[0]
    astorefluc[isn,m]=a
    a=polyfit(log(seps[imin:imax]),log(full[imin:imax,m]),1)[0]
    astorefull[isn,m]=a
  isn=isn+1
  #figure(100)
  #loglog(lam,full[:,3])

avsf=exp(avsf)
#avsf=avsf/nsftot
#avsf=avsf/nsn
avfull=sum(avsf[:,mode,:,:]*nsftot[:,mode,:,:],axis=1)/sum(1.0*nsftot[:,mode,:,:],axis=1)
avperp=avsf[:,mode,-1,:]
avfluc=avsf[:,mode,0,:]
aperps=zeros(nm)
aflucs=zeros(nm)
afulls=zeros(nm)
afull=zeros(nm)
siperps=zeros(nm)
siflucs=zeros(nm)
sifulls=zeros(nm)
for m in xrange(0,nm):
  aperps[m]=polyfit(log(lam[iminperp:imax]),log(avperp[iminperp:imax,m]),1)[0]
  aflucs[m]=polyfit(log(seps[iminperp:imax]),log(avfluc[iminperp:imax,m]),1)[0]
  afulls[m]=polyfit(log(seps[imin:imax]),log(avfull[imin:imax,m]),1)[0]

aperps=mean(astoreperp,axis=0)
aflucs=mean(astorefluc,axis=0)
afulls=mean(astorefull,axis=0)
errfact=1./sqrt(nsn)
for m in xrange(nm):
  siperps[m]=errfact*sqrt(sum((aperps[m]-astoreperp[:,m])**2)/(1.0*nsn))
  siflucs[m]=errfact*sqrt(sum((aflucs[m]-astorefluc[:,m])**2)/(1.0*nsn))
  sifulls[m]=errfact*sqrt(sum((afulls[m]-astorefull[:,m])**2)/(1.0*nsn))

#errperp=array([mean(astoreperp[:5,:],axis=0),mean(astoreperp[5:,:],axis=0)])
#errfluc=array([mean(astorefluc[:5,:],axis=0),mean(astorefluc[5:,:],axis=0)])
#errfull=array([mean(astorefull[:5,:],axis=0),mean(astorefull[5:,:],axis=0)])

#errperptop=zeros(nm)
#errperpbot=zeros(nm)
#errfluctop=zeros(nm)
#errflucbot=zeros(nm)
#errfulltop=zeros(nm)
#errfullbot=zeros(nm)

#for m in xrange(nm):
#  errperptop[m]=max(errperp[:,m])-aperps[m]
#  errperpbot[m]=aperps[m]-min(errperp[:,m])
#  errfluctop[m]=max(errfluc[:,m])-aflucs[m]
#  errflucbot[m]=aflucs[m]-min(errfluc[:,m])
#  errfulltop[m]=max(errfull[:,m])-afulls[m]
#  errfullbot[m]=afulls[m]-min(errfull[:,m])

figure()
loglog(lam,avperp)
for m in xrange(nm):
  loglog(lam,(avperp[imax,m]/lam[imax]**aperps[m])*lam**aperps[m],'k:')
figure()
loglog(lam,avfluc)
for m in xrange(nm):
  loglog(lam,(avfluc[imax,m]/lam[imax]**aflucs[m])*lam**aflucs[m],'k:')
figure()
loglog(lam,avfull)
for m in xrange(nm):
    loglog(lam,(avfull[imax,m]/lam[imax]**afulls[m])*lam**afulls[m],'k:')
mstore=arange(0.0,nm/2.+.5,.5)
mstore2=arange(0.0,10.0,.5)
figure()
plot(mstore2, 1-0.691**mstore2,'k:')
#plot(mstore2,1-sqrt(0.5)**mstore2,'k:')
plot(mstore2,mstore2*(1-sqrt(0.5)**mstore2)/(mstore2/2+1-sqrt(0.5)**mstore2),'k:')
#errorbar(mstore,concatenate([[0],aperps]),yerr=concatenate([[0],siperps]),color='b',linewidth=1.5)
#errorbar(mstore,concatenate([[0],aflucs]),yerr=concatenate([[0],siflucs]),color='g',linewidth=1.5)
#errorbar(mstore,concatenate([[0],afulls]),yerr=concatenate([[0],sifulls]),color='k',linewidth=1.5)

errorbar(mstore,concatenate([[0],aperps]),yerr=[concatenate([[0],siperps]),concatenate([[0],siperps])],color='b',linewidth=1.5)
#errorbar(mstore,concatenate([[0],aflucs]),yerr=[concatenate([[0],siflucs]),concatenate([[0],siflucs])],color='g',linewidth=1.5)
errorbar(mstore,concatenate([[0],afulls]),yerr=[concatenate([[0],sifulls]),concatenate([[0],sifulls])],color='k',linewidth=1.5)
errorbar(ddd[:,0],ddd[:,1],yerr=ddd[:,2],color='r')
xlim(0,10.5)
ylim(0,1)
text(0.5,1.75,r'Total',color='k',fontsize=20)
text(0.5,1.6,r'$\perp$',color='b',fontsize=20)
text(0.5,1.45,r'Fluctuation',color='g',fontsize=20)
title(str(imin)+','+str(imax))
show()
print afulls-ddd[:,1]
