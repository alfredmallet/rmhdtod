aperps=si[-1,-1]
aflucs=si[-1,0]
apars=si[0,-1]
ll=array([0.01,10])
figure(figsize=array([9.5,8]))
loglog(lam,sf[:,0,-1],'r',linewidth=1.5,label=r'$0\degree\leq\theta_{B_{loc}}\!<10\degree,\,80\degree\leq\theta_{\delta z^+}\!<90\degree$')
loglog(lam,sf[:,-1,0],'g',linewidth=1.5,label=r'$80\degree\leq\theta_{B_{loc}}\!<90\degree,\,0\degree\leq\theta_{\delta z^+}\!<10\degree$')
loglog(lam,sf[:,-1,-1],'b',linewidth=1.5,label=r'$80\degree\leq\theta_{B_{loc}}\!<90\degree,\,80\degree\leq\theta_{\delta z^+}\!<90\degree$')
loglog(lam,sf[:,-1,1],'c:',label=r'$80\degree\leq\theta_{B_{loc}}\!<90\degree$')
loglog(lam,sf[:,-1,2:8],'c:')
loglog(lam,sf[:,1,0],'m:',label=r'$0\degree\leq\theta_{\delta z^+}\!<10\degree$')
loglog(lam,sf[:,2:8,0],'m:')

loglog([lam[10],lam[10]],[0.1,10000],'k--')
loglog([lam[27],lam[27]],[0.1,10000],'k--')

loglog(ll,18*ll**aperps,'k')
loglog(ll,13*ll**apars,'k')
loglog(ll,15*ll**aflucs,'k')

aperpstr=str(round(aperps,2))
aflucstr=str(round(aflucs,2))
aparstr=str(round(apars,2))
text(0.25,10,r'$r^{'+aperpstr+'}$',fontsize=20)
text(0.025,0.8,r'$r^{'+aflucstr+'}$',fontsize=20)
text(0.25,2.2,r'$r^{'+aparstr+'}$',fontsize=20)

xlabel(r'$r$',fontsize=30)
ylabel(r'$\langle (\delta z^+)^2 | \theta_{B_{loc}},\theta_{\delta z^+}\rangle$',fontsize=30)

#legend(bbox_to_anchor=(-0.07,0.15,0.9,0.2),frameon=0,fontsize=20)
legend(loc='best',frameon=0,fontsize=20)
xlim(0.02,2)
ylim(0.1,20)
text(1,0.15,'(A)',fontsize=35)
show()
savefig('/Users/alfy/3danis_2014/s2allang.eps')
