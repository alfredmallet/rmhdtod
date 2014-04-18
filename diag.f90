module diag

implicit none

contains

subroutine savesnap(filename,zp,zm)
    use init, only: nlx,nly_par,nlz_par
    use mp, only: iproc
    implicit none
    integer :: i,j,k
    character(len=100),intent(in) :: filename
    real, dimension(nlx,nly_par,nlz_par),intent(in) :: zp,zm
    character(len=10) :: procstr
    write(procstr,"(I0)") iproc
    open(20,file="data/"//trim(procstr)//"/"//trim(filename))
    do k=1,nlz_par
        do j=1,nly_par
            do i=1,nlx
                write(20,*) zp(i,j,k), zm(i,j,k)
            end do
        end do
    end do
    close(20)
end subroutine savesnap

subroutine saveksnap(filename,zpk,zmk)
    use init, only: nky,nkx_par,nlz_par
    use mp, only: iproc
    implicit none
    integer :: i,j,k
    character(len=100),intent(in) :: filename
    complex, dimension(nky,nkx_par,nlz_par),intent(in) :: zpk,zmk
    character(len=10) :: procstr
    write(procstr,"(I0)") iproc
    open(20,file="data/"//trim(procstr)//"/"//trim(filename))
    do k=1,nlz_par
        do j=1,nky
            do i=1,nkx_par
                write(20,99) real(zpk(j,i,k)),aimag(zpk(j,i,k)),real(zmk(j,i,k)),aimag(zmk(j,i,k))
            end do
        end do
    end do
99    format(4g16.8)
    close(20)
end subroutine saveksnap


end module diag
