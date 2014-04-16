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

  nn_to=0
  nn_from=0
  do k=1,nlz
    do j=1,nly
      do i=1,nlx
        if (idx_local(r_variable,j,k)) &
          nn_from(proc_id(k_variable,i,k))=nn_from(proc_id(k_variable,i,k))+1
        if (idx_local(k_variable,i,k)) &
          nn_to(proc_id(r_variable,j,k))=nn_to(proc_id(r_variable,j,k)+1


end module transforms
