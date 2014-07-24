program vecwrite

!*************** Calculates vector fields from scalar fields *************************

use init
use mp
use transforms
use grid
use diag
use subs
implicit none

integer :: i,j,k
character(len=100) :: scalzpf="snap0004.dat",scalzp="/work/01022/tg802750/tod/newtests/force_test/data",vfilename="vec.dat",procstr,runname,inputfile
real, dimension(:,:,:), allocatable :: zp,zm
real, dimension(:,:,:,:), allocatable :: gzp,gzm,vzp,vzm
complex, dimension(:,:,:), allocatable :: zpk,zmk

call getarg(1,runname)
inputfile=trim(runname)//".in"
call read_parameters(inputfile)

!allocate stuff

allocate(zp(nlx,nly_par,nlz_par))
allocate(zm(nlx,nly_par,nlz_par))
allocate(gzp(nlx,nly_par,nlz_par,2))
allocate(gzm(nlx,nly_par,nlz_par,2))
allocate(vzp(nlx,nly_par,nlz_par,2))
allocate(vzm(nlx,nly_par,nlz_par,2))

allocate(zpk(nky,nkx_par,nlz_par))
allocate(zmk(nky,nkx_par,nlz_par))

call init_mp
call init_grid
call init_transforms

call loadsnap(scalzp,scalzpf,zp,zm)

!for testing: load in sine wave

!do k=1,nlz_par
!  do j=1,nly_par
!    do i=1,nlx
!      zp(i,j,k)=ampzp*sin(2.*pi*(kipx*xx(i)/lx+kipy*yy(j)/ly+kipz*zz(k)/lz))
!      zm(i,j,k)=ampzm*sin(2.*pi*(kimx*xx(i)/lx+kimy*yy(j)/ly+kimz*zz(k)/lz))
!    enddo
!  enddo
!enddo

call barrier
!take gradients
do k=1,nlz_par
    call fft(zp(:,:,k),zpk(:,:,k))
    call fft(zm(:,:,k),zmk(:,:,k))
enddo
call grad(zpk,gzp)
call grad(zmk,gzm)

!get vz
vzp(:,:,:,1)=-gzp(:,:,:,2)
vzp(:,:,:,2)=gzp(:,:,:,1)
vzm(:,:,:,1)=-gzm(:,:,:,2)
vzm(:,:,:,2)=gzm(:,:,:,1)
!write vecs to file
write(procstr,"(I0)") iproc
open(10,file=trim(datadir)//"/"//trim(procstr)//"/"//trim(vfilename))
do k=1,nlz_par
    do j=1,nly_par
        do i=1,nlx
            write(*,*) vzp(i,j,k,1), vzp(i,j,k,2), vzm(i,j,k,1), vzm(i,j,k,2)
            write(10,'(4G25.8)') vzp(i,j,k,1), vzp(i,j,k,2), vzm(i,j,k,1), vzm(i,j,k,2)
        enddo
    enddo
enddo
close(10)

call finish_mp

end program vecwrite
