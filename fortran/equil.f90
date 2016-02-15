module equil

!**** Subroutines to do with equilibria ****

use init

implicit none

contains

subroutine makeeq(equilfunc,aeq,deq,sp,sm)

    use init
    use mp, only: iproc, proc0
    implicit none

    integer :: i,j
    character(len=100),intent(in) :: equilfunc
    real, dimension(nlx,nly_par) :: sp,sm
    real :: aeq,deq

    select case (equilfunc)
        case("tear")
            do j=1,nly_par
                do i=1,nlx
                    sp(i,j)=aeq/cosh(xx(i)/deq)**2
                    sm(i,j)=-aeq/cosh(xx(i)/deq)**2
                enddo
            enddo
        case("kh")
            do j=1,nly_par
                do i=1,nlx
                    sp(i,j)=aeq/cosh(xx(i)/deq)**2
                    sm(i,j)=aeq/cosh(xx(i)/deq)**2
                enddo
            enddo

end module equil
