module grid
!***** Contains information about the grids in real and fourier space *****

interface idx_local
  module procedure idx_local_r, idx_local_k
end interface

interface proc_id
  module procedure proc_id_r,proc_id_k
end interface

type :: r_layout_type
   integer :: iproc
   integer :: nlx, nly, nly_par
   integer :: llim_world, ulim_world, llim_proc, ulim_proc
end type r_layout_type

type :: k_layout_type
   integer :: iproc
   integer :: nkx, nky, NKx_par
   integer :: llim_world, ulim_world, llim_proc, ulim_proc
end type k_layout_type

type (r_layout_type) :: r_variable
type (k_layout_type) :: k_variable

contains

subroutine init_grid

  use mp,only: iproc
  use init
  implicit none
  logical, save :: initialized = .false.

  r_variable%iproc = iproc
  r_variable%nlx = nlx
  r_variable%nly = nly
  r_variable%nly_par = nly_par
  r_variable%llim_world = 1
  r_variable%ulim_world = NLy
  r_variable%llim_proc = 1
  r_variable%ulim_proc = NLy_par

  k_variable%iproc = iproc
  k_variable%nkx = nkx
  k_variable%nky = nky
  k_variable%nkx_par = nkx_par
  k_variable%llim_world = 1
  k_variable%ulim_world = NKx
  k_variable%llim_proc = 1
  k_variable%ulim_proc = NKx_par

end subroutine init_grid

subroutine test_grid
  use init
  use mp, only: iproc
  implicit none
  integer :: i,j,k
  if (iproc==0) then
  open(30,file="rgrids.dat")
  open(40,file="kgrids.dat")
  do k=1,nlz
    do j=1,nly
      do i=1,nkx
        write(30,*) i,j,k,proc_id(r_variable,j,k)
        write(40,*) i,j,k,proc_id(k_variable,i,k)
      end do
    end do
  end do
  close(30)
  close(40)
  endif
end subroutine test_grid

function proc_id_r(r_variable,j,zk)
! returns proc no. which has j,zk for a realspace variable
  use init, only:nlz_par, npperp
  implicit none
  integer :: proc_id_r
  type (r_layout_type), intent(in) :: r_variable
  integer, intent(in) :: j,zk
  proc_id_r = (j-1)/r_variable%nly_par + (zk-1)/nlz_par*npperp
end function proc_id_r

function proc_id_k(k_variable,i,zk)
! returns proc no. which has i,zk for a kspace variable
  use init, only: nlz_par,npperp
  implicit none
  integer :: proc_id_k
  type (k_layout_type), intent(in) :: k_variable
  integer, intent(in) :: i,zk
  proc_id_k = (i-1)/k_variable%nkx_par + (zk-1)/nlz_par*npperp
end function proc_id_k

function idx_local_r(r,j,zk)
! returns bool based on whether the proc contains j,zk for realspace variable
  implicit none
  logical :: idx_local_r
  type (r_layout_type), intent(in) :: r
  integer, intent(in) :: j,zk
  idx_local_r = r%iproc == proc_id(r,j,zk)
end function idx_local_r

function idx_local_k (k, i, zk)
! returns bool based on whether the proc contains i,zk for kspace variable
  implicit none
  logical :: idx_local_k
  type (k_layout_type), intent (in) :: k
  integer, intent (in) :: i, zk
  idx_local_k = k%iproc == proc_id(k, i, zk)
end function idx_local_k

real function xx(i)
! returns x value corresponding to i index in the real space grid
    use init
    implicit none
    integer :: i
    xx=lx*(i-1.-nlx/2)/nlx
end function xx

real function yy(j)
! returns y value corresponding to j index in the real space grid
    use init
    use mp, only: iproc
    implicit none
    integer :: j
    yy=ly*(jglobal(j)-1.-nly/2)/nly
end function yy

real function zz(k)
! returns z value corresponding to k index in the real space grid
    use init
    use mp, only: iproc
    implicit none
    integer :: k
    zz=lz*(kglobal(k)-1.-nlz/2)/nlz
end function zz

integer function jglobal(j)
    use init, only: nly_par,npperp
    use mp, only: iproc
    implicit none
    integer :: j
    jglobal=j+mod(iproc,npperp)*nly_par
end function jglobal

integer function kglobal(k)
    use init, only: nlz_par,npperp,npz
    use mp, only: iproc
    implicit none
    integer :: k
    kglobal=k+(iproc/npperp)*nlz_par
end function kglobal

real function kx(i)
  use init, only: npperp,nkx_par
  use mp, only: iproc
  implicit none
  integer :: i,iglobal
  iglobal = mod(iproc,npperp)*nkx_par+i
  kx=iglobal-1.
end function kx

real function ky(j)
  use init, only: nky,lx,ly
  implicit none
  integer :: j
  if (j<=nky/2+1) then
    ky=(j-1.)*lx/ly
  else
    ky=(j-nky-1.)*lx/ly
  endif
end function ky

real function kperp(j,i)
  implicit none
  integer :: j,i
  kperp=sqrt(ky(j)**2+kx(i)**2)
end function kperp

end module grid
