module transforms

use init
use fft_work, only: fft_type, POINTER_KIND
use redistribute, only: redist_type

implicit none

type (redist_type), save :: r2k
type (fft_type) :: fft_x2k, fft_k2x, fft_y2k, fft_k2y!

integer(Kind=POINTER_KIND) :: iplan_x2k, iplan_k2x, iplan_y2k, iplan_k2y

contains

subroutine init_transforms

  implicit none

  call init_redistribute

end subroutine init_transforms

subroutine init_redistribute

    use grid
    use mp
    use redistribute, only: index_list_type,init_redist,delete_list
    implicit none
    
    type (index_list_type), dimension(0:nproc-1) :: to_list, from_list
    integer, dimension(0:nproc-1) :: nn_to, nn_from
    integer, dimension (2) :: from_low, from_high, to_high
    integer :: to_low
    logical :: initialized = .false.
    integer :: i, j, ip, n, k
  
    if (initialized) return
    initialized=.true.
    
    !counts elements that need to go to/from each proc
    nn_to=0
    nn_from=0
    do k=1,nlz
        do j=1,nly
            do i=1,nkx
                if (idx_local(r_variable,j,k)) &
                nn_from(proc_id(k_variable,i,k))=nn_from(proc_id(k_variable,i,k))+1
                if (idx_local(k_variable,i,k)) &
                nn_to(proc_id(r_variable,j,k))=nn_to(proc_id(r_variable,j,k))+1
            enddo
        enddo
    enddo
  
    do ip=0,nproc-1
        if (nn_from(ip)>0) then
            allocate(from_list(ip)%first(nn_from(ip)))
            allocate(from_list(ip)%second(nn_from(ip)))
        endif
        if (nn_to(ip)>0) then
            allocate(to_list(ip)%first(nn_to(ip)))
            allocate(to_list(ip)%second(nn_to(ip)))
        endif
    enddo
  
    !get local indices of the elements that will go to/from other procs
    nn_to=0
    nn_from=0
    do k=1, nlz
        do j = 1, nly
            do i = 1, nkx   ! this is dealiasing in x
                if (idx_local(r_variable,j,k)) then
                    ip = proc_id(k_variable,i,k)
                    n = nn_from(ip)+1
                    nn_from(ip) = n
                    from_list(ip)%first(n) = i
                    from_list(ip)%second(n) = 1+mod(j-1,nly_par)
                end if
                if (idx_local(k_variable,i,k)) then
                    ip = proc_id(r_variable,j,k)
                    n = nn_to(ip)+1
                    nn_to(ip) = n
                    to_list(ip)%first(n) = j
                    to_list(ip)%second(n) = 1+mod(i-1,nkx_par)
                end if
            enddo
        enddo
    enddo
  
    from_low(1)=1
    from_low(2)=1
    to_low=1
    to_high(1)=nly
    to_high(2)=nkx_par
    from_high(1)=nlx/2+1
    from_high(2)=nly_par
  
    call init_redist(r2k,'c',to_low,to_high,to_list,from_low,from_high,from_list)

    call delete_list(to_list)
    call delete_list(from_list)

end subroutine init_redistribute

!********** Subroutines for doing forward transforms **********

subroutine fft(array,arrayk)
! 2d transform from real to k space
    use fft_work
    use redistribute, only: gather
    use mp
    use grid, only: kx,ky,kperp
    implicit none
  
    real, dimension(:,:) :: array
    complex, dimension(:,:) :: arrayk
    complex, allocatable, dimension(:,:) :: array_temp, ak
    integer :: i,j
  
    allocate(array_temp(nlx/2+1,nly_par))
    allocate(ak(nly,nkx_par))
  
    array_temp=0.0
    ak=0.0
  
    call xfft(array,array_temp)
    ! Hou-Li filtering
    do i=1,nkx
        array_temp(i,:)=array_temp(i,:)*exp(-36.0*((i*1.0-1.0)/((nkx-1)*1.0))**36)
    enddo
  
    call gather(r2k,array_temp,ak)

    call yfft(ak)
    ! Hou-Li filtering
    do j=1,nky
        arrayk(j,:)=ak(j,:)*exp(-36.0*(abs(ky(j))/(nky/2.*lx/ly))**36)
    enddo

    deallocate(ak)
    deallocate(array_temp)

end subroutine fft

subroutine xfft(array,ak)
! does the fft in the x direction
    use fft_work
    implicit none
  
    real, dimension(:,:) :: array
    complex, dimension(:,:) :: ak
    logical, save :: first = .true.
    real, save :: scale
    integer :: imax, jmax
  
    imax=size(array,1)
    jmax=size(array,2)

    if (first) then
        call init_rcfftw(fft_x2k,-1,imax,iplan_x2k)
        scale=1./sqrt(real(imax))
        first = .false.
    endif

    call rfftwnd_f77_real_to_complex(iplan_x2k,jmax,array,1,imax,ak,1,imax/2+1)
    ak=ak*scale

end subroutine xfft

subroutine yfft(ak)
! does the fft in the y direction
    use fft_work
    implicit none
    
    complex, dimension(:,:) :: ak
    logical, save :: first = .true.
    complex :: dummy
    real, save :: scale
    integer :: imax,jmax

    imax=size(ak,1)
    jmax=size(ak,2)

    if (first) then
        iplan_y2k=10 ! do fft in place
        call init_ccfftw(fft_y2k,-1,imax,iplan_y2k)
        scale = 1./sqrt(real(imax))
        first=.false.
    end if

    call fftwnd_f77(iplan_y2k,jmax,ak,1,imax,dummy,0,0)
    ak=ak*scale

end subroutine yfft

!********** Subroutines for doing inverse transforms **********

subroutine ifft(arrayk,array)
! does the inverse 2d transform    
    use fft_work
    use mp, only: iproc
    use redistribute, only: scatter
    implicit none

    real, dimension(:,:) :: array
    complex, dimension(:,:) :: arrayk
    complex, allocatable, dimension(:,:) :: array_temp,ak
    integer :: i,j,iglobal

    allocate (ak(nly,nkx_par))
    ak=arrayk
    call iyfft(ak)

    allocate(array_temp(nlx/2+1,nly_par))
    array_temp=0.0
    call scatter(r2k,ak,array_temp)
    call ixfft(array_temp,array)

    deallocate(array_temp)
    deallocate(ak)

end subroutine ifft

subroutine ixfft(akx,array)
! inverse x transform    
    use fft_work
    implicit none

    real, dimension(:,:) :: array
    complex, dimension(:,:) :: akx
    logical, save :: first=.true.
    real, save :: scale
    integer :: i,j
    
    i=size(array,1)
    j=size(array,2)
    
    if (first) then
        call init_crfftw(fft_k2x,1,i,iplan_k2x)
        scale=1./sqrt(real(i))
        first=.false.
    endif

    call rfftwnd_f77_complex_to_real(iplan_k2x,j,akx,1,i/2+1,array,1,i)
    array=array*scale

end subroutine ixfft

subroutine iyfft(ak)
! inverse y transform    
    use fft_work
    implicit none

    complex, dimension(:,:) :: ak
    logical, save :: first=.true.
    complex :: dummy
    real, save :: scale
    integer :: i,j,k1,k2

    i=size(ak,1)
    j=size(ak,2)

    if (first) then
        iplan_k2y=10 !do fft in place
        call init_ccfftw(fft_k2y,1,i,iplan_k2y)
        scale=1./sqrt(real(i))
        first=.false.
    endif

    call fftwnd_f77(iplan_k2y,j,ak,1,i,dummy,0,0)
    ak=ak*scale

end subroutine iyfft

end module transforms
