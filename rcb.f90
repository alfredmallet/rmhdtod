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
integer :: isep,isamp,nsep=32,nsamp=10000,i0,j0,k0,skp,ifile
real :: phi,theta,di,dj,dk,ti,tj,tk,azp,azm,au,ab,thetaBloc,thetapm,thetaub
real, dimension(3) :: Bloc,Blochat
real, dimension(2) :: dzp,dzm,du,db,zpp,zpm,zmp,zmm,dr
real, dimension(2) :: zpmpar,zmmpar,zpppar,zmppar,dzmh_par,dzph_par
real, dimension(:), allocatable :: seps
logical :: zp_found,zm_found
integer :: diz_int,i1p_0,i1p_1,j1p_0,j1p_1,i1m_0,i1m_1,j1m_0,j1m_1,km,kp
integer :: im_0,im_1,jm_0,jm_1,ip_0,ip_1,jp_0,jp_1
real :: diypar,dixpar,dzmh,dzph,dzph_old,dzmh_old,lzp,lparp,lzm,lparm,tmx,tmy,tpx,tpy
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

do ifile=0,npperp*orig_npz-1
    
    write(*,*) "Processor", iproc, "reporting that I am reading data in", ifile

    write(fold0,"(I0)") ifile
    open(10,file=trim(vfold)//"/"//trim(fold0)//"/"//trim(vfilename))
    do k=1,zpointsperfile
        do j=1,nly_par
            do i=1,nlx
                read(10,'(4G25.8)') vzp(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,1), vzp(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,2), vzm(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,1), vzm(i,j+mod(ifile,npperp)*nly_par,k+(ifile/npperp)*zpointsperfile,2)
            enddo
        enddo
    enddo
    close(10)
enddo
write(*,*) iproc, 'done reading'

!correct bug
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
        !Separation
        di=cos(phi)*seps(isep)
        dj=sin(phi)*seps(isep)
        dr=(/di,dj/)
        !Find final points
        i1p_0=modulo(i0+floor(di/2.),nlx)
        i1p_1=modulo(i0+ceiling(di/2.),nlx)
        j1p_0=modulo(j0+floor(dj/2.),nly)
        j1p_1=modulo(j0+ceiling(dj/2.),nly)
        if (i1p_0.eq.0) i1p_0=nlx
        if (i1p_1.eq.0) i1p_1=nlx
        if (j1p_0.eq.0) j1p_0=nly
        if (j1p_1.eq.0) j1p_1=nly
        i1m_0=modulo(i0-floor(di/2.),nlx)
        i1m_1=modulo(i0-ceiling(di/2.),nlx)
        j1m_0=modulo(j0-floor(dj/2.),nly)
        j1m_1=modulo(j0-ceiling(dj/2.),nly)
        if (i1m_0.eq.0) i1m_0=nlx
        if (i1m_1.eq.0) i1m_1=nlx
        if (j1m_0.eq.0) j1m_0=nly
        if (j1m_1.eq.0) j1m_1=nly
        !Fraction along
        ti=di-floor(di)
        tj=dj-floor(dj)
        zpp=(1-ti)*(1-tj)*vzp(i1p_0,j1p_0,k0,:)+ti*(1-tj)*vzp(i1p_1,j1p_0,k0,:)&
           +(1-ti)*(tj)*vzp(i1p_0,j1p_1,k0,:)+(ti)*(tj)*vzp(i1p_1,j1p_1,k0,:)
        zmp=(1-ti)*(1-tj)*vzm(i1p_0,j1p_0,k0,:)+ti*(1-tj)*vzm(i1p_1,j1p_0,k0,:)&
           +(1-ti)*(tj)*vzm(i1p_0,j1p_1,k0,:)+(ti)*(tj)*vzm(i1p_1,j1p_1,k0,:)
        zpm=(1-ti)*(1-tj)*vzp(i1m_0,j1m_0,k0,:)+ti*(1-tj)*vzp(i1m_1,j1m_0,k0,:)&
           +(1-ti)*(tj)*vzp(i1m_0,j1m_1,k0,:)+(ti)*(tj)*vzp(i1m_1,j1m_1,k0,:)
        zmm=(1-ti)*(1-tj)*vzm(i1m_0,j1m_0,k0,:)+ti*(1-tj)*vzm(i1m_1,j1m_0,k0,:)&
           +(1-ti)*(tj)*vzm(i1m_0,j1m_1,k0,:)+(ti)*(tj)*vzm(i1m_1,j1m_1,k0,:)
        
        !Calculate deltas
        dzp(:)=zpp-zpm
        dzm(:)=zmp-zmm
        du=0.5*(dzp+dzm)
        db=0.5*(dzp-dzm)
        !calculate amplitudes
        azp=sqrt(dzp(1)**2+dzp(2)**2)
        azm=sqrt(dzm(1)**2+dzm(2)**2)
        au=sqrt(du(1)**2+du(2)**2)
        ab=sqrt(db(1)**2+db(2)**2)
        !calculate direction of local mean field
        Bloc=(/0d0,0d0,1d0/)+0.25*(/zpp(1)-zmp(1)+zpm(1)-zmm(1),zpp(2)-zmp(2)+zpm(2)-zmm(2),0.0/)
        Blochat=Bloc/sqrt(Bloc(1)**2+Bloc(2)**2+Bloc(3)**2)
        !calculate angle between separation and local mean field
        thetaBloc=acos((dr(1)*blochat(1)+dr(2)*blochat(2)+dr(3)*blochat(3))/sqrt(dr(1)**2+dr(2)**2+dr(3)**2))
        !calculate angles between flucts
        thetapm=asin((dzp(1)*dzm(2)-dzp(2)*dzm(1))/sqrt((dzp(1)**2+dzp(2)**2)*(dzm(1)**2+dzm(2)**2)))
        thetaub=asin((du(1)*db(2)-du(2)*db(1))/sqrt((du(1)**2+du(2)**2)*(db(1)**2+db(2)**2)))
        
        !finding parallel correlation length
        dzph=0.
        dzmh=0.
        dixpar=0.
        diypar=0.
        diz_int=0.
        zp_found=.false.
        zm_found=.false.
        do
          if (.not.zp_found) then
            lzp=nlz*1.0
            lparp=nlz*1.0/abs(Blochat(3))
          endif
          if (.not.zm_found) then
            lzm=nlz*1.0
            lparm=nlz*1.0/abs(Blochat(3))
          endif
          diz_int=diz_int+2
          dixpar=dixpar+2.0*Blochat(1)/Blochat(3)
          diypar=diypar+2.0*Blochat(2)/Blochat(3)
          km=modulo(k0-diz_int/2,nlz)
          kp=modulo(k0+diz_int/2,nlz)
          if (km.eq.0) km=nlz
          if (kp.eq.0) kp=nlz
          im_0=modulo(i0-floor(dixpar/2.),nlx)
          im_1=modulo(i0-ceiling(dixpar/2.),nlx)
          jm_0=modulo(j0-floor(diypar/2.),nly)
          jm_1=modulo(j0-ceiling(diypar/2.),nly)
          tmx=dixpar/2.-floor(dixpar/2.)
          tmy=diypar/2.-floor(diypar/2.)
          if (im_0.eq.0) im_0=nlx
          if (im_1.eq.0) im_1=nlx
          if (jm_0.eq.0) jm_0=nly
          if (jm_1.eq.0) jm_1=nly
          ip_0=modulo(i0+floor(dixpar/2.),nlx)
          ip_1=modulo(i0+ceiling(dixpar/2.),nlx)
          jp_0=modulo(j0+floor(diypar/2.),nly)
          jp_1=modulo(j0+ceiling(diypar/2.),nly)
          tpx=dixpar/2.-floor(dixpar/2.)
          tpy=diypar/2.-floor(diypar/2.)
          if (ip_0.eq.0) ip_0=nlx
          if (ip_1.eq.0) ip_1=nlx
          if (jp_0.eq.0) jp_0=nly
          if (jp_1.eq.0) jp_1=nly
          zpmpar=(1-tmx)*((1-tmy)*vzp(im_0,jm_0,km,:)+tmy*vzp(im_0,jm_1,km,:)) &
                 +tmx*((1-tmy)*vzp(im_1,jm_0,km,:)+tmy*vzp(im_1,jm_1,km,:))
          zmmpar=(1-tmx)*((1-tmy)*vzm(im_0,jm_0,km,:)+tmy*vzm(im_0,jm_1,km,:)) &
                 +tmx*((1-tmy)*vzm(im_1,jm_0,km,:)+tmy*vzm(im_1,jm_1,km,:))
          zpppar=(1-tpx)*((1-tpy)*vzp(ip_0,jp_0,kp,:)+tpy*vzp(ip_0,jp_1,kp,:)) &
                 +tpx*((1-tpy)*vzp(ip_1,jp_0,kp,:)+tpy*vzp(ip_1,jp_1,kp,:))
          zmppar=(1-tpx)*((1-tpy)*vzp(ip_0,jp_0,kp,:)+tpy*vzp(ip_0,jp_1,kp,:)) &
                 +tpx*((1-tpy)*vzp(ip_1,jp_0,kp,:)+tpy*vzp(ip_1,jp_1,kp,:))
          dzph_old=dzph
          dzmh_old=dzmh
          dzph_par=zpmpar-zpppar
          dzmh_par=zmmpar-zmppar
          dzph=sqrt(dzph_par(1)**2+dzph_par(2)**2)
          dzmh=sqrt(dzmh_par(1)**2+dzmh_par(2)**2)
          if ((dzph>azp).and.(.not.zp_found)) then
            zp_found=.true.
            lzp=diz_int*1.0-2.+2.*((azp-dzph_old)/(dzph-dzph_old))
            lparp=lzp/abs(Blochat(3))
          endif
          if ((dzmh>azm).and.(.not.zm_found)) then
            zm_found=.true.
            lzm=diz_int*1.0-2.+2.*((azm-dzmh_old)/(dzmh-dzmh_old))
            lparm=lzm/abs(Blochat(3))
          endif
          if ((zp_found.and.zm_found).or.(diz_int.ge.nlz)) exit
        enddo
        !write everything to files
        write(35,"(11G25.8)") azp,azm,au,ab,thetaBloc,thetapm,thetaub,lparp,lparm,lzp,lzm
    enddo
    close(35)
    enddo

call finish_mp
end program sf
