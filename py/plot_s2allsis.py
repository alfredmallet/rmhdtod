from os import listdir

wdir=0
wsf=0
m=3
iminpar=10
imaxpar=27
iminperp=10
imaxperp=30
lsd=listdir('/Users/alfy/1024_sfdata')
nsn=1.*len(lsd)
seps=zeros(32)
isn=0
sibit=zeros((9,9,int16(nsn)))
sf=zeros((32,9,9))
for isep in xrange(32):
    seps[isep]=4.*(1024.*0.25/4)**((isep)*1./31.)
lam=seps*pi/512.

for lss in lsd:
    print isn,lss
    d=loadtxt('/Users/alfy/1024_sfdata/'+lss)
    sfbit=d[:,6+wsf]
    sfbit=sfbit.reshape((32,9,4,9,10))[:,:,wdir,:,m]
    sf=sf+log(sfbit)/nsn
    for ang in xrange(9):
        for ali in xrange(9):
            if ang==0:
                imin=iminpar
                imax=imaxpar
            else:
                imin=iminperp
                imax=imaxperp
            sibit[ang,ali,isn]=polyfit(log(lam[imin:imax]),log(sfbit[imin:imax,ang,ali]),1)[0]
    isn=isn+1

sf=exp(sf)
si=zeros((9,9))
err=zeros((9,9))


for ang in xrange(9):
    for ali in xrange(9):
        if ang==0:
            imin=iminpar
            imax=imaxpar
        else:
            imin=iminperp
            imax=imaxperp
        si[ang,ali]=polyfit(log(lam[imin:imax]),log(sf[imin:imax,ang,ali]),1)[0]
        err[ang,ali]=sqrt(sum((si[ang,ali]-sibit[ang,ali,:])**2))*2./nsn

figure(figsize=array([9.5,8]))
imshow(si,interpolation='nearest',origin='lower')
xlim(-0.5,8.4999)
for ang in xrange(9):
    for ali in xrange(9):
        text(-0.31+ali,-0.27+ang,r'$'+str('%.2f' % si[ang,ali])+'$\n'+r'$\!\pm'+str((round(err[ang,ali],2)))[1:]+'$',color='white',fontsize=15)

xticks([0,2,4,6,8],[5,25,45,65,85])
yticks([0,2,4,6,8],[5,25,45,65,85])

if m%2:
    mstr=str(m/2.+0.5)
else:
    mstr=m/2
sfstr=['(\delta z^+_\perp)^'+mstr,'(\delta z^-_\perp)^'+mstr,'\delta b_\perp^'+mstr,'\delta u_\perp^'+mstr]
dirstr=['\delta z^+_\perp','\delta z^-_\perp','\delta b_\perp','\delta u_\perp']

xlabel(r'$\theta_{'+dirstr[wdir]+'}$',fontsize=30)
ylabel(r'$\theta_{B_{loc}}$',fontsize=30)
text(-1.5,-1.1,'(B)',fontsize=35)
savefig('s2allsis.eps')
