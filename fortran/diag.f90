module diag

implicit none

integer :: ispecfile,isnapfile

contains

subroutine savesnap(filename,zp,zm,t)

    use init, only: nlx,nly_par,nlz_par,datadir
    use mp, only: iproc,proc0
    implicit none
    
    integer :: i,j,k
    character(len=100),intent(in) :: filename
    real, dimension(nlx,nly_par,nlz_par),intent(in) :: zp,zm
    character(len=10) :: procstr
    real :: t

    if (proc0) then
        open(25,file=trim(datadir)//"/"//"tsnaps.dat",access="append")
        write(25,*) isnapfile,t
        close(25)
    endif
    write(procstr,"(I0)") iproc
    open(20,file=trim(datadir)//"/"//trim(procstr)//"/"//trim(filename))
    do k=1,nlz_par
        do j=1,nly_par
            do i=1,nlx
                write(20,*) zp(i,j,k), zm(i,j,k)
            enddo
        enddo
    enddo
    close(20)

end subroutine savesnap

subroutine loadsnap(locsnap,fsnap,zp,zm)
    use init, only: nlx,nly_par,nlz_par
    use mp, only: iproc,proc0
    implicit none

    integer :: i,j,k
    character(len=100),intent(in) :: locsnap,fsnap
    real, dimension(nlx,nly_par,nlz_par) :: zp,zm
    character(len=10) :: procstr
    
    write(procstr,"(I0)") iproc
    open(75,file=trim(locsnap)//"/"//trim(procstr)//"/"//trim(fsnap))
    do k=1,nlz_par
        do j=1,nly_par
            do i=1,nlx
                read(75,*) zp(i,j,k), zm(i,j,k)
            enddo
        enddo
    enddo
    close(75)

end subroutine loadsnap

subroutine loadeq(loceq,feq,zpeq,zmeq)
    use init, only: nlx,nly_par,nlz_par
    use mp, only: iproc,proc0
    implicit none
    
    integer :: i,j,k
    character(len=100),intent(in) :: loceq,feq
    real,dimension(nlx,nly_par) :: zpeq,zmeq
    character(len=10) :: procstr

    write(procstr,"(I0)") iproc
    open(76,file=trim(loceq)//"/"//trim(procstr)//"/"//trim(feq))
    do j=1,nly_par
        do i=1,nlx
            read(76,*) zpeq(i,j), zmeq(i,j)
        enddo
    enddo
    close(76)

end subroutine loadeq

subroutine saveksnap(filename,zpk,zmk)
!   not used: may be useful for tests    
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
            enddo
        enddo
    enddo
99  format(4g16.8)
    close(20)

end subroutine saveksnap

end module diag
