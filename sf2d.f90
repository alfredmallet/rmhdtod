program sf

!****** Calculates structure function increments *******
! Should be run on 32 large mem cpus
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
integer :: isep,isamp,nsep=32,nsamp=1000,i0,j0,k0,i1_0,i1_1,j1_0,j1_1,f1_000,f1_010,f1_001,f1_011,j1_0_loc,j1_1_loc,skp,ifile,izfile,ipfile
real :: phi,theta,di,dj,dk,ti,tj,tk,zpx1,zpy1,zmx1,zmy1,avzpx1,avzpy1,avzmx1,avzmy1,azp,azm,au,ab,thetaBloc,rdotBloc,zpdotBloc,zmdotBloc,azpperp,azmperp,auperp,abperp,thetazp,thetazm,thetau,thetab,thetapm,thetaub,thetapmperp,thetaubperp
real, dimension(2) :: rperp,rhatperp,zpperp,zmperp,uperp,bperp,zpperphat,zmperphat,uperphat,bperphat
real, dimension(2) :: dr,dzp,dzm,du,db
real, dimension(:), allocatable :: seps

integer :: cc,cm
real :: cr

call getarg(1,runname)
inputfile=trim(runname)//".in"

call getarg(2,sfdir)
call getarg(3,vfold)

call read_parameters(inputfile)

orig_npz=nlz/zpointsperfile

allocate(vzp(nlx,nly,zpointsperfile,2))
allocate(vzm(nlx,nly,zpointsperfile,2))
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

do izfile=0,orig_npz-1
do ipfile=0,npperp-1
    ifile=ipfile+izfile*npperp
    write(*,*) "Processor", iproc, "reporting that I am reading data in", ifile

    write(fold0,"(I0)") ifile
    open(10,file=trim(vfold)//"/"//trim(fold0)//"/"//trim(vfilename))
    do k=1,zpointsperfile
        do j=1,nly_par
            do i=1,nlx
                read(10,*) vzp(i,j+mod(ifile,npperp)*nly_par,k,1), vzp(i,j+mod(ifile,npperp)*nly_par,k,2), vzm(i,j+mod(ifile,npperp)*nly_par,k,1), vzm(i,j+mod(ifile,npperp)*nly_par,k,2)
            enddo
        enddo
    enddo
    close(10)
enddo
write(*,*) iproc, 'done reading for plane ', izfile
!Loop over separations
do isep=1,nsep
    write(*,*) iproc,isep,izfile
    write(wfile,"(I0)") isep
    write(wfold,"(I0)") iproc
    open(35,file=trim(sfdir)//"/"//trim(wfold)//"/base"//trim(wfile)//".dat",access='append')
    !Loop to get random samples
    do isamp=1,nsamp
        !Random initial point
        i0=nlx*.9999*uniran()+1
        j0=nly*.9999*uniran()+1
        k0=zpointsperfile*.9999*uniran()+1
        !Random angles
        phi=pi*(uniran()*2.0)
        !Separation
        di=cos(phi)*seps(isep)
        dj=sin(phi)*seps(isep)
        dr=(/di,dj/)
        !Find final points
        i1_0=modulo(i0+floor(di),nlx)
        i1_1=modulo(i0+ceiling(di),nlx)
        j1_0=modulo(j0+floor(dj),nly)
        j1_1=modulo(j0+ceiling(dj),nly)
        if (i1_0.eq.0) i1_0=nlx
        if (i1_1.eq.0) i1_1=nlx
        if (j1_0.eq.0) j1_0=nly
        if (j1_1.eq.0) j1_1=nly
        !Fraction along
        ti=di-floor(di)
        tj=dj-floor(dj)
        !000
        avzpx1=vzp(i1_0,j1_0,k0,1)*(1-ti)*(1-tj)
        avzpy1=vzp(i1_0,j1_0,k0,2)*(1-ti)*(1-tj)
        avzmx1=vzm(i1_0,j1_0,k0,1)*(1-ti)*(1-tj)
        avzmy1=vzm(i1_0,j1_0,k0,2)*(1-ti)*(1-tj)
        !100
        avzpx1=avzpx1+vzp(i1_1,j1_0,k0,1)*ti*(1-tj)
        avzpy1=avzpy1+vzp(i1_0,j1_0,k0,2)*ti*(1-tj)
        avzmx1=avzmx1+vzm(i1_1,j1_0,k0,1)*ti*(1-tj)
        avzmy1=avzmy1+vzm(i1_1,j1_0,k0,2)*ti*(1-tj)
        !010
        avzpx1=avzpx1+vzp(i1_0,j1_1,k0,1)*(1-ti)*(tj)
        avzpy1=avzpy1+vzp(i1_0,j1_1,k0,2)*(1-ti)*(tj)
        avzmx1=avzmx1+vzm(i1_0,j1_1,k0,1)*(1-ti)*(tj)
        avzmy1=avzmy1+vzm(i1_0,j1_1,k0,2)*(1-ti)*(tj)
        !110
        avzpx1=avzpx1+vzp(i1_1,j1_1,k0,1)*(ti)*(tj)
        avzpy1=avzpy1+vzp(i1_1,j1_1,k0,2)*(ti)*(tj)
        avzmx1=avzmx1+vzm(i1_1,j1_1,k0,1)*(ti)*(tj)
        avzmy1=avzmy1+vzm(i1_1,j1_1,k0,2)*(ti)*(tj)
       
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
        rperp=dr
        rhatperp=rperp/sqrt(rperp(1)**2+rperp(2)**2+rperp(3)**2)
        !calculate local perpendicular fluctuation directions
        zpperp=dzp
        zmperp=dzm
        uperp=0.5*(zpperp+zmperp)        
        bperp=0.5*(zpperp-zmperp)
        zpperphat=zpperp/sqrt(zpperp(1)**2+zpperp(2)**2+zpperp(3)**2)
        zmperphat=zmperp/sqrt(zmperp(1)**2+zmperp(2)**2+zmperp(3)**2)
        uperphat=uperp/sqrt(uperp(1)**2+uperp(2)**2+uperp(3)**2)
        bperphat=bperp/sqrt(bperp(1)**2+bperp(2)**2+bperp(3)**2)

        !calculate angle between perp separation and perp local fluctuations
        thetazp=acos(rhatperp(1)*zpperphat(1)+rhatperp(2)*zpperphat(2))
        thetazm=acos(rhatperp(1)*zmperphat(1)+rhatperp(2)*zmperphat(2))
        thetau=acos(rhatperp(1)*uperphat(1)+rhatperp(2)*uperphat(2))
        thetab=acos(rhatperp(1)*bperphat(1)+rhatperp(2)*bperphat(2))
        !calculate angles between flucts
        thetapm=acos((dzp(1)*dzm(1)+dzp(2)*dzm(2))/sqrt((dzp(1)**2+dzp(2)**2)*(dzm(1)**2+dzm(2)**2)))
        thetaub=acos((du(1)*db(1)+du(2)*db(2))/sqrt((du(1)**2+du(2)**2)*(db(1)**2+db(2)**2)))
        !write everything to files
        write(35,"(10G25.8)") azp,azm,au,ab,thetazp,thetazm,thetau,thetab,thetapm,thetaub
    enddo
    close(35)
    enddo
enddo
call finish_mp
end program sf
