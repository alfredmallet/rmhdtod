from rgbl import *
rgbl=get_rgbl(nm)
figure(figsize=array([9.5,8]))
alamlpar=zeros(nm)
lsd=listdir('1024_sfdata')
nsn=len(lsd)
mode=0
mode2=0
iminperp=10
imaxperp=27
imin=iminperp
imax=imaxperp
isn=0
ast=zeros((nsn,nm))
si=zeros(nm)
perpav=zeros((32,nm))
flucav=zeros((32,nm))
parav=zeros((32,nm))
for lss in lsd:
    print isn, lss
    d=loadtxt(ddir+'/'+lss)
    zpsf=d[:,6+mode]
    print zpsf[0]
    nsf=d[:,5]
    nsfr=nsf.reshape((32,9,4,9,nm))
    zpsfr=zpsf.reshape((32,9,4,9,nm))
    perp=zpsfr[:,-1,mode2,-1,:]
    fluc=zpsfr[:,-1,mode2,0,:]
    par=zpsfr[:,0,mode2,-1,:]
    perpav=perpav+log(perp)/nsn
    flucav=flucav+log(fluc)/nsn
    parav=parav+log(par)/nsn
    for m in xrange(nm):
      ast[isn,m]=polyfit(log(par[imin:imax,m]),log(perp[imin:imax,m]),1)[0]
    isn=isn+1
avperp=exp(perpav)
avfluc=exp(flucav)
avpar=exp(parav)
for i in xrange(nm):
    alamlpar[i]=polyfit(log(avpar[iminperp:imaxperp,i]),log(avperp[iminperp:imaxperp,i]),1)[0]
    loglog(avpar[:,i],avperp[:,i],color=rgbl[nm-1-i],linewidth=1.5)
    errfact=2./sqrt(nsn*1.0)
    si[i]=errfact*sqrt(sum((alamlpar[i]-ast[:,i])**2)/(1.0*nsn))
#loglog(avpar[iminperp,:],avperp[iminperp,:],'k')
#loglog(avpar[imaxperp,:],avperp[imaxperp,:],'k')

loglog([1e-15*ones(nm),avpar[15,:],1e15*ones(nm)],avperp[15,:]*(array([1e-15*ones(nm),avpar[15,:],1e15*ones(nm)])/avpar[15,:])**alamlpar,'k:')
#loglog([1e-15*ones(nm),avpar[15,:],1e15*ones(nm)],avperp[15,:]*(array([1e-15*ones(nm),avpar[15,:],1e15*ones(nm)])/avpar[15,:])**0.5,'k:')

xlabel(r'$\langle (\delta z^+_\perp)^n | \theta_{B_{loc}} = 0\degree \rangle$',fontsize=30)
ylabel(r'$\langle (\delta z^+_\perp)^n | \theta_{B_{loc}} = 90\degree,\theta_{\delta z^+} = 90\degree \rangle$',fontsize=30)
text(150,3e3,r'$n=5.0$',color=rgbl[0],fontsize=30)
text(1.0,0.6,r'$n=0.5$',color=rgbl[-1],fontsize=30)
xlim(0.1,1e4)
ylim(0.1,1e4)
ax=gca()
ax.tick_params(axis='both',which='both',labelsize=20)
tight_layout()

a=axes([0.59,0.22,0.33,0.33])
errorbar(ms,1-alamlpar,yerr=si,color='k',linewidth=1.5)
xlabel(r'$n$',fontsize=25)
ylabel(r'$\gamma_n$',fontsize=25)
ylim(0.2,0.5)

savefig('/Users/alfy/3danis_2014/lamlpar.eps')
