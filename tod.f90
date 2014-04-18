program tod
!********************************************************************************
!*                          TOD - solves RMHD equations                         *
!*                  praise be to the ancestors, GOSTA and VIRIATO               *
!********************************************************************************

use init
use mp
use transforms
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

!********** Initialization **********

call getarg(1,runname)
inputfile=trim(runname)//".in"
call read_parameters(inputfile)

allocate(zp(nlx,nly_par,nlz_par))
allocate(zm(nlx,nly_par,nlz_par))

allocate(zpk(nky,nkx_par,nlz_par))
allocate(zmk(nky,nkx_par,nlz_par))

call init_mp
call init_grid
call init_transforms

!***** Setting initial fields, if required *****

if (initfield.eq."wave") then
    if (proc0) then
        write(*,*) "Sinusoidal initial field with given wavenumbers and amplitudes, kip=",kipx,kipy,kipz,"kim=",kimx,kimy,kimz,"ampzp=",ampzp,"ampzm=",ampzm
    endif
    do k=1,nlz_par
        do j=1,nly_par
            do i=1,nlx
                zp(i,j,k)=ampzp*sin(2.*pi*(kipx*xx(i)/lx+kipy*yy(j)/ly+kipz*zz(k)/lz))
                zm(i,j,k)=ampzm*sin(2.*pi*(kimx*xx(i)/lx+kimy*yy(j)/ly+kimz*zz(k)/lz))
            end do
        end do
    end do
else
    zp(:,:,:)=0
    zm(:,:,:)=0
endif

!TODO: other initfield options, like "norm"
!TEMP: testing forward fft
do k=1,nlz_par
  call fft(zp(:,:,k),zpk(:,:,k))
  call fft(zm(:,:,k),zmk(:,:,k))
end do
do ip=0,nproc-1
  if (iproc.eq.ip) then
    write(*,*) iproc,size(zpk,1),size(zpk,2),size(zpk,3)
  endif
  call barrier
end do

!start of the timestep
it=0

if ((ifields.gt.0).and.(mod(it,ifields).eq.0)) then
    write(itstr,"(I0)") it
    filename="snap"//trim(itstr)//".dat"
    call savesnap(filename,zp,zm)
    call saveksnap("k"//filename,zpk,zmk)
endif

end program tod
