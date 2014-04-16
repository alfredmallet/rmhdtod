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

end module diag
