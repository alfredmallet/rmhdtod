module grid
!***** Contains information about the grids in real and fourier space *****

!temp until mpi is supported

contains

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
