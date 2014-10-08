rgbl=get_rgbl(nm)
figure(figsize=array([9.5,8]))
axilpar=zeros(nm)
for i in xrange(nm):
    iminperp=10
    imaxperp=30
    axilpar[i]=polyfit(log(avpar[iminperp:imaxperp,i]),log(avfluc[iminperp:imaxperp,i]),1)[0]
    loglog(avpar[:,i],avfluc[:,i],color=rgbl[nm-1-i],linewidth=1.5)

loglog(avpar[iminperp,:],avfluc[iminperp,:],'k')
loglog(avpar[imaxperp,:],avfluc[imaxperp,:],'k')

loglog([1e-15*ones(nm),avpar[15,:],1e15*ones(nm)],avfluc[15,:]*(array([1e-15*ones(nm),avpar[15,:],1e15*ones(nm)])/avpar[15,:])**axilpar,'k:')
#loglog([1e-15*ones(nm),avpar[15,:],1e15*ones(nm)],avfluc[15,:]*(array([1e-15*ones(nm),avpar[15,:],1e15*ones(nm)])/avpar[15,:])**0.5,'k:')

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
plot(ms,1-axilpar,'k')
xlabel(r'$n$',fontsize=25)
ylabel(r'$\gamma_n$',fontsize=25)
ylim(0.3,1)

savefig('/Users/alfy/3danis_2014/xilpar.eps')
