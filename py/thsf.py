#run 1024chiplotter first
figure(figsize=array([9.5,8]))
rgbl=get_rgbl(10)
for i in xrange(10):
    loglog(lam,thsf[:,i],color=rgbl[9-i],linewidth=1.5)
    loglog(lam,thsf[-1,i]*(lam/lam[-1])**((1-alamxi[i])*ms[i]),'k')

xlim(0.02,5)
xlabel(r'$r$',fontsize=30)
ylabel(r'$\sin^n\!{\theta}^{ub}_{n}$',fontsize=30)
text(2.0,0.5,r'$n=0.5$',color=rgbl[-1],fontsize=30)
text(1.5,0.02,r'$n=5.0$',color=rgbl[0],fontsize=30)
ax=gca()
a=axes([0.55,0.15,0.33,0.33])
xlabel(r'$n$',fontsize=25)
ylabel(r'$n\alpha_n$',fontsize=25)
plot(ms,athsf,'r')
plot(ms,(1-alamxi)*ms,'k')
