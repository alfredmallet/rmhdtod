module init
!***** Contains all information required to initialize the code *****
implicit none
save

!***** Declare constants *****

real, parameter :: pi=3.1415926535
!box
real :: lx=1,ly=1,lz=1
integer :: nlx=16,nly=16,nlz=16
integer :: npperp=1,npz=1
!time
real :: tmax=1
integer :: imax=100
real :: cfl_frac=0.25
!diss
integer :: hyper_order=3
real :: nu=0.0
!start
real :: ampzp=1.0,ampzm=0.5
character(len=4) :: initfield="norm"
real, dimension(3) :: kip=(/1.0,0.0,1.0/),kim=(/0.0,1.0,1.0/)
!force
logical :: turb=.true.
real :: kfp1=1,kfp2=2,kfz1=1,kfz2=1,epsm=0.0,epsp=-1
!output
integer :: iout=50,ispec=-1,ifields=1000000
real :: tspec=2.*pi,tfields=2.*pi
!restart
integer :: restart=0
character(len=100) :: rsfile="old_fields.dat"
character(len=100) :: rspath=""
!definitions
integer :: nkx, nky
integer :: kmax,kperpmax
integer :: nkx_par,nly_par,nlz_par
real :: dx,dy,dz
real :: scale

!note on physics: Alfv\'en speed v_A is degenerate with lz in the RMHD equations.
!To change v_A, change lz.
contains

subroutine read_parameters(inputfile)
    !Reads parameters from the input file
    implicit none
    character(len=100), intent(in) :: inputfile
    integer  :: ierr

    namelist /box_parameters/ lx,ly,lz,nlx,nly,nlz,npperp,npz
    namelist /time_parameters/ tmax, imax,cfl_frac
    namelist /diss_parameters/ hyper_order,nu
    namelist /start_parameters/ ampzp,ampzm,initfield,kip,kim
    namelist /force_parameters/ turb,kfp1,kfp2,kfz1,kfz2,epsp,epsm
    namelist /output_parameters/ iout, ispec, ifields, tspec, tfields
    namelist /restart_parameters/ restart,rsfile,rspath

    open(unit=10,file=trim(inputfile),status='old')
    read(10,nml=box_parameters,iostat=ierr)
        if(ierr /= 0) then
            write(*,*) "Reading box_parameters failed"
            write(*,*) "lx,ly,lz,nlx,nly,nlz,npperp,npz=", lx,ly,lz,nlx,nly,nlz,npperp,npz
        endif
    read(10,nml=time_parameters,iostat=ierr)
        if(ierr /= 0) then
            write(*,*) "Reading time_parameters failed"
            write(*,*) "tmax, imax,cfl_frac=", tmax, imax,cfl_frac
        endif
    read(10,nml=diss_parameters,iostat=ierr)
        if(ierr /= 0) then
            write(*,*) "Reading diss_parameters failed"
            write(*,*) "hyper_order,nu=", hyper_order,nu
        endif
    read(10,nml=start_parameters,iostat=ierr)
        !if(ierr /= 0) then
         !   write(*,*) "Reading start_parameters failed"
          !  write(*,*) "ampzp,ampzm,initfield,kip,kim=", ampzp,ampzm,initfield,kip,kim
        !endif
    read(10,nml=force_parameters,iostat=ierr)
        if(ierr /= 0) then
            write(*,*) "Reading force_parameters failed"
            write(*,*) "turb,kfp1,kfp2,kfz1,kfz2,epsp,epsm=", turb,kfp1,kfp2,kfz1,kfz2,epsp,epsm
        endif
    read(10,nml=output_parameters,iostat=ierr)
        if(ierr /= 0) then
            write(*,*) "Reading output_parameters failed"
            write(*,*) "iout, ispec, ifields, tspec, tfields=", iout, ispec, ifields, tspec, tfields
        endif
    read(10,nml=restart_parameters,iostat=ierr)
        if (ierr/=0) then
            write(*,*) "Reading restart_parameters failed"
            write(*,*) "restart,rsfile,rspath", restart,rsfile,rspath
        endif
    close(10)

    lx=2.*pi*lx
    ly=2.*pi*ly
    lz=2.*pi*lz

    !***** Definitions *****
    nkx=nlx/2+1
    nky=nly
    kmax=nkx-nkx/2 !maximum k after filtering
    kperpmax=sqrt((nkx*1.0)**2+(nky/2.)**2)
    nkx_par=(nkx-1)/npperp+1
    nly_par=(nly-1)/npperp+1
    nlz_par=(nlz-1)/npz+1
    dx=lx/(nlx*1.)
    dy=ly/(nly*1.)
    dz=lz/(nlz*1.)
    scale=1./(nlx*nly) ! for ffts
end subroutine read_parameters

end module init
