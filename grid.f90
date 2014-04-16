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
   integer :: NLx, NLy, NLy_par
   integer :: llim_world, ulim_world, llim_proc, ulim_proc
end type r_layout_type

type :: k_layout_type
   integer :: iproc
   integer :: NKx, NKy, NKx_par
   integer :: llim_world, ulim_world, llim_proc, ulim_proc
end type k_layout_type

type (r_layout_type) :: r_variable
type (k_layout_type) :: k_variable

contains

subroutine init_grid

  use mp,only: iproc
  use constants
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

function proc_id_r(r_variable,j,zk)
! returns proc no. which has j,zk for a realspace variable
  use constants, only:nlz_par, NPE
  implicit none
  integer :: proc_id_r
  type (r_layout_type),intent(in) :: r_variable

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

end module grid
