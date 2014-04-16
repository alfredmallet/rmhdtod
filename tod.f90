program tod
!********************************************************************************
!*                          TOD - solves RMHD equations                         *
!*                  praise be to the ancestors, GOSTA and VIRIATO               *
!********************************************************************************

use init
use mp
!use transforms
use grid
use diag
!use forcing

implicit none
integer :: i,j,k,it,ip

real, dimension(:,:,:), allocatable :: zp,zm
complex,dimension(:,:,:), allocatable :: zpk,zmk

character(len=100) :: runname, inputfile
character(len=100) :: filename
character(len=10) :: itstr
character(len=10) :: procstr
character(len=100) :: datadir="data"

!temporary - until mpi is done
!write(procstr,"(I0)") iproc
!write(*,*) "procstr",procstr
!********** Initialization **********
call init_mp
call getarg(1,runname)
inputfile=trim(runname)//".in"
call read_parameters(inputfile)
if (proc0) then
    write(*,*) "Read in parameters: lx,ly,lz,nlx,nly,nlz,npperp,npz=",lx,ly,lz,nlx,nly,nlz,npperp,npz
endif
allocate(zp(nlx,nly_par,nlz_par))
allocate(zm(nlx,nly_par,nlz_par))

allocate(zpk(nky,nkx_par,nlz_par))
allocate(zmk(nky,nkx_par,nlz_par))
write(*,*) nlx,nly_par,nlz_par
if (initfield.eq."wave") then
    if (proc0) then
        write(*,*) "Sinusoidal initial field with given wavenumbers and amplitudes"
    endif
    do k=1,nlz_par
        do j=1,nly_par
            do i=1,nlx
                zp(i,j,k)=ampzp*sin(2*pi*kip(3)*zz(k)/lz)
!(2.*pi*(kip(1)*xx(i)/lx+kip(2)*yy(j)/ly+kip(3)*zz(k)/lz))
                zm(i,j,k)=ampzm*sin(2.*pi*(kim(1)*xx(i)/lx+kim(2)*yy(j)/ly+kim(3)*zz(k)/lz))
            end do
        end do
    end do
else
    zp(:,:,:)=0
    zm(:,:,:)=0
endif

!TODO: other initfield options, like "norm"
call barrier

do ip=0,npperp*npz-1
  if (iproc==ip) then
    write(*,*) "iproc",ip,"xx"
    do i=1,nlx
    write(*,*) xx(i)/lx
    end do
  endif
  call barrier
end do

do ip=0,npperp*npz-1
  if (iproc==ip) then
    write(*,*) "iproc",ip,"yy"
    do j=1,nly_par
    write(*,*) yy(j)/ly
    end do
  endif
  call barrier
end do
do ip=0,npperp*npz-1
  if (iproc==ip) then
    write(*,*) "iproc",ip, "zz"
    do k=1,nlz_par
    write(*,*) zz(k)/lz,sin(2*pi*zz(k)/lz),zp(1,1,k)
    end do
  endif
  call barrier
end do

!start of the timestep
it=0

if ((ifields.gt.0).and.(mod(it,ifields).eq.0)) then
    write(itstr,"(I0)") it
    filename="snap"//trim(itstr)//".dat"
    call savesnap(filename,zp,zm)
endif

end program tod
