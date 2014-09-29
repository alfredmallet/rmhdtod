program sf

!****** Calculates structure function increments *******
! Should be run on ~16 large mem cpus for 1024^3
! Change npz to 1 in input file, set zpointsperfile here
use init
use mp
use transforms
use grid
use diag
use subs
implicit none

integer :: orig_npz,zpointsperfile=4
integer :: i,j,k
character (len=100) :: fold0,vfilename="vec.dat",sfdir,wfold,wfile,runname,inputfile,vfold
real, dimension(:,:,:,:), allocatable :: vzp,vzm
integer :: isep,isamp,nsep=32,nsamp=500000,i0,j0,k0,i1_0,i1_1,j1_0,j1_1,k1_0,k1_1,f1_000,f1_010,f1_001,f1_011,j1_0_loc,j1_1_loc,k1_0_loc,k1_1_loc,skp,ifile
real :: phi,theta,di,dj,dk,ti,tj,tk,zpx1,zpy1,zmx1,zmy1,avzpx1,avzpy1,avzmx1,avzmy1,azp,azm,au,ab,thetaBloc,rdotBloc,zpdotBloc,zmdotBloc,azpperp,azmperp,auperp,abperp,thetazp,thetazm,thetau,thetab,thetapm,thetaub,thetapmperp,thetaubperp
real, dimension(3) :: dr,Bloc,Blochat,rperp,rhatperp,zpperp,zmperp,uperp,bperp,zpperphat,zmperphat,uperphat,bperphat
real, dimension(2) :: dzp,dzm,du,db
real, dimension(:), allocatable :: seps

integer :: cc,cm
real :: cr

call getarg(1,runname)
inputfile=trim(runname)//".in"
call read_parameters(inputfile)

call getarg(2,sfdir)
call getarg(3,vfold)

orig_npz=nlz/zpointsperfile

allocate(vzp(nlx,nly,nlz,2))
allocate(vzm(nlx,nly,nlz,2))
allocate(seps(nsep))

call init_mp
call system_clock(cc,cr,cm)
cr=uniran(init=cc+iproc*floor(2e8))
call barrier
do isep=1,nsep
  seps(isep)=4.0*(nlx*0.25/4.0)**((isep-1)*1.0/(nsep-1))
  if (iproc==0) write(*,*) seps(isep)
  write(wfile,"(I0)") isep
  write(wfold,"(I0)") iproc
  open(35,file=trim(sfdir)//"/"//trim(wfold)//"/base"//trim(wfile)//".dat")
  close(35)
enddo

write(*,*) "NLX,NLY,NLZ=", nlx,nly,nlz

do ifile=0,npperp*orig_npz-1
    
    write(*,*) "Processor", iproc, "reporting that I am reading data in", ifile

    write(fold0,"(I0)") ifile
    open(10,file=trim(vfold)//"/"//trim(fold0)//"/"//trim(vfilename))
    do k=1,zpointsperfile
        do j=1,nly_par
            do i=1,nlx
                read(10,'(4G25.8)') vzp(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,1), vzp(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,2), vzm(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,1), vzm(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,2)
                !if ((ifile>1022).and.(iproc.eq.0)) write(*,*) i, j+mod(ifile,npperp)*nly_par, k+(ifile/npperp)*zpointsperfile, vzp(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,:),vzm(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,:)
            enddo
        enddo
    enddo
    close(10)
enddo
write(*,*) iproc, 'done reading'

!correct bug (I HOPE...)
vzp=-vzp
vzm=-vzm

!Loop over separations
do isep=1,nsep
    write(*,*) iproc,isep
    write(wfile,"(I0)") isep
    write(wfold,"(I0)") iproc
    open(35,file=trim(sfdir)//"/"//trim(wfold)//"/base"//trim(wfile)//".dat",access='append')
    !Loop to get random samples
    do isamp=1,nsamp
        !Random initial point
        i0=nlx*.9999*uniran()+1
        j0=nly*.9999*uniran()+1
        k0=nlz*.9999*uniran()+1
        !Random angles
        phi=pi*(uniran()*2.0)
        theta=pi*(uniran()-0.5)
        !Separation
        di=sin(theta)*cos(phi)*seps(isep)
        dj=sin(theta)*sin(phi)*seps(isep)
        dk=cos(theta)*seps(isep)
        dr=(/di,dj,dk/)
        !Find final points
        i1_0=modulo(i0+floor(di),nlx)
        i1_1=modulo(i0+ceiling(di),nlx)
        j1_0=modulo(j0+floor(dj),nly)
        j1_1=modulo(j0+ceiling(dj),nly)
        k1_0=modulo(k0+floor(dk),nlz)
        k1_1=modulo(k0+ceiling(dk),nlz)
        if (i1_0.eq.0) i1_0=nlx
        if (i1_1.eq.0) i1_1=nlx
        if (j1_0.eq.0) j1_0=nly
        if (j1_1.eq.0) j1_1=nly
        if (k1_0.eq.0) k1_0=nlz
        if (k1_1.eq.0) k1_1=nlz
        !Fraction along
        ti=di-floor(di)
        tj=dj-floor(dj)
        tk=dk-floor(dk)
        !000
        avzpx1=vzp(i1_0,j1_0,k1_0,1)*(1-ti)*(1-tj)*(1-tk)
        avzpy1=vzp(i1_0,j1_0,k1_0,2)*(1-ti)*(1-tj)*(1-tk)
        avzmx1=vzm(i1_0,j1_0,k1_0,1)*(1-ti)*(1-tj)*(1-tk)
        avzmy1=vzm(i1_0,j1_0,k1_0,2)*(1-ti)*(1-tj)*(1-tk)
        !100
        avzpx1=avzpx1+vzp(i1_1,j1_0,k1_0,1)*ti*(1-tj)*(1-tk)
        avzpy1=avzpy1+vzp(i1_0,j1_0,k1_0,2)*ti*(1-tj)*(1-tk)
        avzmx1=avzmx1+vzm(i1_1,j1_0,k1_0,1)*ti*(1-tj)*(1-tk)
        avzmy1=avzmy1+vzm(i1_1,j1_0,k1_0,2)*ti*(1-tj)*(1-tk)
        !010
        avzpx1=avzpx1+vzp(i1_0,j1_1,k1_0,1)*(1-ti)*(tj)*(1-tk)
        avzpy1=avzpy1+vzp(i1_0,j1_1,k1_0,2)*(1-ti)*(tj)*(1-tk)
        avzmx1=avzmx1+vzm(i1_0,j1_1,k1_0,1)*(1-ti)*(tj)*(1-tk)
        avzmy1=avzmy1+vzm(i1_0,j1_1,k1_0,2)*(1-ti)*(tj)*(1-tk)
        !110
        avzpx1=avzpx1+vzp(i1_1,j1_1,k1_0,1)*(ti)*(tj)*(1-tk)
        avzpy1=avzpy1+vzp(i1_1,j1_1,k1_0,2)*(ti)*(tj)*(1-tk)
        avzmx1=avzmx1+vzm(i1_1,j1_1,k1_0,1)*(ti)*(tj)*(1-tk)
        avzmy1=avzmy1+vzm(i1_1,j1_1,k1_0,2)*(ti)*(tj)*(1-tk)
        !001
        avzpx1=avzpx1+vzp(i1_0,j1_0,k1_1,1)*(1-ti)*(1-tj)*(tk)
        avzpy1=avzpy1+vzp(i1_0,j1_0,k1_1,2)*(1-ti)*(1-tj)*(tk)
        avzmx1=avzmx1+vzm(i1_0,j1_0,k1_1,1)*(1-ti)*(1-tj)*(tk)
        avzmy1=avzmy1+vzm(i1_0,j1_0,k1_1,2)*(1-ti)*(1-tj)*(tk)
        !101
        avzpx1=avzpx1+vzp(i1_1,j1_0,k1_1,1)*(ti)*(1-tj)*(tk)
        avzpy1=avzpy1+vzp(i1_1,j1_0,k1_1,2)*(ti)*(1-tj)*(tk)
        avzmx1=avzmx1+vzm(i1_1,j1_0,k1_1,1)*(ti)*(1-tj)*(tk)
        avzmy1=avzmy1+vzm(i1_1,j1_0,k1_1,2)*(ti)*(1-tj)*(tk) 
        !011
        avzpx1=avzpx1+vzp(i1_0,j1_1,k1_1,1)*(1-ti)*(tj)*(tk)
        avzpy1=avzpy1+vzp(i1_0,j1_1,k1_1,2)*(1-ti)*(tj)*(tk)
        avzmx1=avzmx1+vzm(i1_0,j1_1,k1_1,1)*(1-ti)*(tj)*(tk)
        avzmy1=avzmy1+vzm(i1_0,j1_1,k1_1,2)*(1-ti)*(tj)*(tk)
        !111
        avzpx1=avzpx1+vzp(i1_1,j1_1,k1_1,1)*(ti)*(tj)*(tk)
        avzpy1=avzpy1+vzp(i1_1,j1_1,k1_1,2)*(ti)*(tj)*(tk)
        avzmx1=avzmx1+vzm(i1_1,j1_1,k1_1,1)*(ti)*(tj)*(tk)
        avzmy1=avzmy1+vzm(i1_1,j1_1,k1_1,2)*(ti)*(tj)*(tk)
        !Calculate deltas
        dzp(:)=vzp(i0,j0,k0,:)-(/avzpx1,avzpy1/)
        dzm(:)=vzm(i0,j0,k0,:)-(/avzmx1,avzmy1/)
        du=0.5*(dzp+dzm)
        db=0.5*(dzp-dzm)
        !calculate amplitudes
        azp=sqrt(dzp(1)**2+dzp(2)**2)
        azm=sqrt(dzm(1)**2+dzm(2)**2)
        au=sqrt(du(1)**2+du(2)**2)
        ab=sqrt(db(1)**2+db(2)**2)
        !write(*,*) iproc,azp,dr,i0,j0,k0,vzp(i0,j0,k0,:),avzpx1,avzpy1,ti,tj,tk
        !calculate direction of local mean field
        Bloc=(/0d0,0d0,1d0/)+0.25*(/(vzp(i0,j0,k0,1)-vzm(i0,j0,k0,1)+avzpx1-avzmx1),(vzp(i0,j0,k0,2)-vzm(i0,j0,k0,2)+avzpy1-avzmy1),0.0/)
        Blochat=Bloc/sqrt(Bloc(1)**2+Bloc(2)**2+Bloc(3)**2)
        !calculate angle between separation and local mean field
        thetaBloc=acos((dr(1)*blochat(1)+dr(2)*blochat(2)+dr(3)*blochat(3))/sqrt(dr(1)**2+dr(2)**2+dr(3)**2))
        !calculate local perp separation
        rdotBloc=dr(1)*blochat(1)+dr(2)*blochat(2)+dr(3)*blochat(3)
        rperp=dr-rdotBloc*Blochat
        rhatperp=rperp/sqrt(rperp(1)**2+rperp(2)**2+rperp(3)**2)
        !calculate local perpendicular fluctuation directions
        zpdotBloc=dzp(1)*blochat(1)+dzp(2)*blochat(2)
        zpperp=(/dzp(1),dzp(2),0.0/)-zpdotBloc*blochat
        zmdotBloc=dzm(1)*blochat(1)+dzm(2)*blochat(2)
        zmperp=(/dzm(1),dzm(2),0.0/)-zmdotBloc*blochat
        uperp=0.5*(zpperp+zmperp)        
        bperp=0.5*(zpperp-zmperp)
        zpperphat=zpperp/sqrt(zpperp(1)**2+zpperp(2)**2+zpperp(3)**2)
        zmperphat=zmperp/sqrt(zmperp(1)**2+zmperp(2)**2+zmperp(3)**2)
        uperphat=uperp/sqrt(uperp(1)**2+uperp(2)**2+uperp(3)**2)
        bperphat=bperp/sqrt(bperp(1)**2+bperp(2)**2+bperp(3)**2)

        !calculate perpendicular amplitudes
        azpperp=sqrt(zpperp(1)**2+zpperp(2)**2+zpperp(3)**2)
        azmperp=sqrt(zmperp(1)**2+zmperp(2)**2+zmperp(3)**2)
        auperp=sqrt(uperp(1)**2+uperp(2)**2+uperp(3)**2)
        abperp=sqrt(bperp(1)**2+bperp(2)**2+bperp(3)**2)
        !calculate angle between perp separation and perp local fluctuations
        thetazp=acos(rhatperp(1)*zpperphat(1)+rhatperp(2)*zpperphat(2)+rhatperp(3)*zpperphat(3))
        thetazm=acos(rhatperp(1)*zmperphat(1)+rhatperp(2)*zmperphat(2)+rhatperp(3)*zmperphat(3))
        thetau=acos(rhatperp(1)*uperphat(1)+rhatperp(2)*uperphat(2)+rhatperp(3)*uperphat(3))
        thetab=acos(rhatperp(1)*bperphat(1)+rhatperp(2)*bperphat(2)+rhatperp(3)*bperphat(3))
        !calculate angles between flucts
        thetapm=acos((dzp(1)*dzm(1)+dzp(2)*dzm(2))/sqrt((dzp(1)**2+dzp(2)**2)*(dzm(1)**2+dzm(2)**2)))
        thetaub=acos((du(1)*db(1)+du(2)*db(2))/sqrt((du(1)**2+du(2)**2)*(db(1)**2+db(2)**2)))
        !calculate perp angles between flucts
        thetapmperp=acos(zpperphat(1)*zmperphat(1)+zpperphat(2)*zmperphat(2)+zpperphat(3)*zmperphat(3))
        thetaubperp=acos(uperphat(1)*bperphat(1)+uperphat(2)*bperphat(2)+uperphat(3)*bperphat(3))

        !write everything to files
        write(35,"(11G25.8)") azp,azm,au,ab,thetaBloc,thetazp,thetazm,thetau,thetab,thetapm,thetaub
    enddo
    close(35)
    enddo

call finish_mp
end program sf
