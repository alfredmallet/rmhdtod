from os import listdir
from numpy import *
import pickle
from mayavi import mlab
from scipy.interpolate import UnivariateSpline as uspl

imin=0
imax=32

sfstr=['z_\perp^+','z_\perp^-','u_\perp','b_\perp']

lsd=listdir('1024_sfdata')
nsn=len(lsd)
minsurf=zeros((4,10,4))
maxsurf=ones((4,10,4))*1000000
pows=zeros((4,10,4,16))
seps=zeros(32)
sf=zeros((32,9,4,9,10,4))
totnsf=zeros((32,9,4,9,10))
r=zeros((16,9*2,9*4,4,10,4))
x=zeros((16,9*2,9*4,4,10,4))
y=zeros((16,9*2,9*4,4,10,4))
z=zeros((16,9*2,9*4,4,10,4))
rr=zeros((16,9*2,9*4+1,4,10,4))
xx=zeros((16,9*2,9*4+1,4,10,4))
yy=zeros((16,9*2,9*4+1,4,10,4))
zz=zeros((16,9*2,9*4+1,4,10,4))
for isep in xrange(32):
    seps[isep]=4.*(1024.*0.25/4)**((isep)*1./31.)
for lss in lsd:
    print lss
    d=loadtxt('1024_sfdata/'+lss)
    zpsf=d[:,6:10]
    nsf=d[:,5]
    nsf=nsf.reshape((32,9,4,9,10))
    zpsf=zpsf.reshape((32,9,4,9,10,4))
    sf=sf+log(zpsf)
    totnsf=totnsf+nsf

sf=exp(sf/len(lsd))

for ang in xrange(9):
    for wdir in xrange(4):
        for ali in xrange(9):
            for m in xrange(10):
                for wsf in xrange(4):
                    minsurf[wdir,m,wsf]=maximum(minsurf[wdir,m,wsf],min(sf[imin:imax,ang,wdir,ali,m,wsf]))
                    maxsurf[wdir,m,wsf]=minimum(maxsurf[wdir,m,wsf],max(sf[imin:imax,ang,wdir,ali,m,wsf]))
minsurf=minsurf
maxsurf=maxsurf
for wdir in xrange(4):
    for m in xrange(10):
        for wsf in xrange(4):
            pows[wdir,m,wsf,:]=minsurf[wdir,m,wsf]*(maxsurf[wdir,m,wsf]/minsurf[wdir,m,wsf])**arange(0.,1+1./15.,1./15.)
for p in xrange(16):
    for ang in xrange(18):
        ang2=ang
        if ang>=9:
            ang2=8-ang%9
        for ali in xrange(36):
            ali2=ali
            if (ali>=9 and ali<18):
                ali2=8-ali%9
            if (ali>=18 and ali<27):
                ali2=ali%9
            if (ali>=27):
                ali2=8-ali%9
            for wdir in xrange(4):
                for m in xrange(10):
                    for wsf in xrange(4):
                        fr=uspl(log(sf[imin:imax,ang2,wdir,ali2,m,wsf]),log(seps[imin:imax]))
                        r[p,ang,ali,wdir,m,wsf]=exp(fr(log(pows[wdir,m,wsf,p])))*pi/512
                        x[p,ang,ali,wdir,m,wsf]=r[p,ang,ali,wdir,m,wsf]*sin(pi*0.5*ang/9.0+0.25/9.0*pi)*cos(pi*0.5*ali/9.0+0.25/9.0*pi)
                        y[p,ang,ali,wdir,m,wsf]=r[p,ang,ali,wdir,m,wsf]*sin(pi*0.5*ang/9.0+0.25/9.0*pi)*sin(pi*0.5*ali/9.0+0.25/9.0*pi)
                        z[p,ang,ali,wdir,m,wsf]=r[p,ang,ali,wdir,m,wsf]*cos(pi*0.5*ang/9.0+0.25/9.0*pi)

rr[:,:,:9*4,:,:,:]=r
rr[:,:,-1,:,:,:]=r[:,:,0,:,:,:]
xx[:,:,:9*4,:,:,:]=x
xx[:,:,-1,:,:,:]=x[:,:,0,:,:,:]
yy[:,:,:9*4,:,:,:]=y
yy[:,:,-1,:,:,:]=y[:,:,0,:,:,:]
zz[:,:,:9*4,:,:,:]=z
zz[:,:,-1,:,:,:]=z[:,:,0,:,:,:]
print rr[0,8,8,0,:,0]

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
for wdir in xrange(1):
    for m in xrange(10):
        for wsf in xrange(1):
            for p in [15]:
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')
                colors=np.empty((9*2,9*4+1,4))
                for ang in xrange(9*2):
                    for ali in xrange(9*4+1):
                        colors[ang,ali,:]=cm.jet((rr[p,ang,ali,wdir,m,wsf]-numpy.min(rr[p,:,:,wdir,m,wsf]))/(numpy.max(rr[p,:,:,wdir,m,wsf])-numpy.min(rr[p,:,:,wdir,m,wsf])))
                ax.plot_surface(xx[p,:,:,wdir,m,wsf], yy[p,:,:,wdir,m,wsf], zz[p,:,:,wdir,m,wsf],  rstride=1, cstride=1, facecolors=colors)
                #xlim(-numpy.max(rr[p,:,:,wdir,m,wsf]),numpy.max(rr[p,:,:,wdir,m,wsf]))
                #ylim(-numpy.max(rr[p,:,:,wdir,m,wsf]),numpy.max(rr[p,:,:,wdir,m,wsf]))
                #ax.set_zlim3d(-numpy.max(rr[p,:,:,wdir,m,wsf]),numpy.max(rr[p,:,:,wdir,m,wsf]))
                xlim(-1.5*rr[p,8,8,wdir,m,wsf],1.5*rr[p,8,8,wdir,m,wsf])
                ylim(-1.5*rr[p,8,8,wdir,m,wsf],1.5*rr[p,8,8,wdir,m,wsf])
                ax.set_zlim3d(-1.5*rr[p,8,8,wdir,m,wsf],1.5*rr[p,8,8,wdir,m,wsf])
                ax.set_aspect('equal')
                ax.azim=45.
                ax.elev=30.
                #plot(xx[p,8,:,wdir,m,wsf],yy[p,8,:,wdir,m,wsf],ones(37)*-numpy.max(rr[p,:,:,wdir,m,wsf]),'k:')
                #plot(ones(37)*-numpy.max(rr[p,:,:,wdir,m,wsf]),concatenate([yy[p,:,8,wdir,m,wsf],-yy[p,:,8,wdir,m,wsf],ones(1)*(yy[p,0,8,wdir,m,wsf])]),concatenate([zz[p,:,8,wdir,m,wsf],-zz[p,:,8,wdir,m,wsf],ones(1)*zz[p,0,8,wdir,m,wsf]]),'k:') 
                #plot(concatenate([xx[p,:,0,wdir,m,wsf],-xx[p,:,0,wdir,m,wsf],ones(1)*(xx[p,0,0,wdir,m,wsf])]),ones(37)*-numpy.max(rr[p,:,:,wdir,m,wsf]),concatenate([zz[p,:,0,wdir,m,wsf],-zz[p,:,0,wdir,m,wsf],ones(1)*zz[p,0,0,wdir,m,wsf]]),'k:')
                plot(xx[p,8,:,wdir,m,wsf],yy[p,8,:,wdir,m,wsf],ones(37)*-1.5*rr[p,8,8,wdir,m,wsf],'k:')
                plot(ones(37)*-1.5*rr[p,8,8,wdir,m,wsf],concatenate([yy[p,:,8,wdir,m,wsf],-yy[p,:,8,wdir,m,wsf],ones(1)*(yy[p,0,8,wdir,m,wsf])]),concatenate([zz[p,:,8,wdir,m,wsf],-zz[p,:,8,wdir,m,wsf],ones(1)*zz[p,0,8,wdir,m,wsf]]),'k:')
                plot(concatenate([xx[p,:,0,wdir,m,wsf],-xx[p,:,0,wdir,m,wsf],ones(1)*(xx[p,0,0,wdir,m,wsf])]),ones(37)*-1.5*rr[p,8,8,wdir,m,wsf],concatenate([zz[p,:,0,wdir,m,wsf],-zz[p,:,0,wdir,m,wsf],ones(1)*zz[p,0,0,wdir,m,wsf]]),'k:')
                tstr=r'$\langle (\delta '+sfstr[wsf]+')^{'+str(m/2.+0.5)+'} \\rangle ='
                title(tstr+'%.2f$' % pows[wdir,m,wsf,p])
                xlabel(r'$\xi$')
                ylabel(r'$\lambda$')
                ax.set_zlabel(r'$l_\parallel$')
                if p<10:
                    pp='0'+str(p)
                else:
                    pp=str(p)
                savefig('1024surfs/'+str(wdir)+str(wsf)+pp+str(m)+'.png')
                close()


f=open('1024surfs/surfdata.dat','w')
pickle.dump([r,pows],f)
f.close()

#doing inferred alignment
lam=r[:,8,8,:,:,:]
xi=r[:,8,0,:,:,:]
lpar=r[:,0,8,:,:,:]
thlamxi=lam/xi
thlaml=lam/lpar
thxil=xi/lpar

lam00=lam[:,0,:,0]
xi00=xi[:,0,:,0]
lpar00=lpar[:,0,:,0]
thlamxi00=thlamxi[:,0,:,0]
thlaml00=thlaml[:,0,:,0]
thxil00=thxil[:,0,:,0]

for i in xrange(10):
    axil[i]=polyfit(log(lam[10:,0,i,0]),log(thxil[10:,0,i,0]),1)[0]
    alamxi[i]=polyfit(log(lam[10:,0,i,0]),log(thlamxi[10:,0,i,0]),1)[0]
    alaml[i]=polyfit(log(lam[10:,0,i,0]),log(thlaml[10:,0,i,0]),1)[0]
