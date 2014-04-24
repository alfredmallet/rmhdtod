module subs

use init

implicit none

contains

subroutine force(kfz,eps,dt,spk,smk,field)

    use init
    use mpi, only: proc0
    implicit none

    real, intent(in) :: dt,eps,kfz
    character, intent(in) :: field
    integer, dimension(:,:), allocatable :: kfs_big
    integer, dimension(:,:), allocatable, save :: kfs
    complex, dimension(:,:,:) :: spk,smk
    real :: kfp, kav=0
    real :: amp, phi, phiz
    integer :: i,j,ik,nk=0,sgn,i_rl,i_im
    logical, save :: lfirst=.true.

    ! calculate possible elements
    if (lfirst) then
        allocate(kfs_big(nkx*nky,2))
        do i=-ceiling(kfp2),ceiling(kfp2)
            do j=-ceiling(kfp2,ceiling(kfp2)
                kfp=sqrt(i**2+j**2)
                    if ((kfp.ge.kfp1).and.(kfp.le.kfp2).and.(kfp.ne.0.)) then
                        nk=nk+1
                        kfs_big(nk,1)=i
                        kfs_big(nk,2)=j
                        kav=kav+kfp
                    endif
            enddo
        enddo
        kav=kav/nk
        allocate(kfs(nk,2))
        kfs=kfs_big(1:nk,1:2)
        deallocate(kfs_big)
        lfirst=.false.
    endif
    
    ! choose random element
    ik=nk*.9999*uniran()+1
    kfp=sqrt(kfs(ik,1)**2.+kfs(ik,2)**2.)
    ! choose amplitude and phase
    amp=sqrt(-eps/dt/kfp**2.0*log(uniran()))*sqrt(nlx*nly) !N.B. may depend on fft convention
    phi=pi*(2.0*uniran()-1.0)
    ! convert to array indices
    sgn=sign(1,kfs(ik,1))
    i=sgn*kfs(ik,1)+1
    j=sgn*kfs(ik,2)*ly/lx+1
    if (j.le.0) j=j+nky
    proc_force=(i-1)/nkx_par
    iloc=i-proc_force*nkx_par
    phi=sgn*phi
    ! stuff for forcing the z-dirn mode
    if (kfz.ne.0) then
        amp=amp*sqrt(2.0)
        phiz=pi*(2.0*uniran()-1.0)
    else
        phiz=0.0
    endif

    if (mod(iproc,npperp).eq.proc_force) then
        if ((field.eq.'b').or.(field.eq.'p') then
            do k=1,nlz_par
                spk(j,iloc,k)=spk(j,iloc,k)+amp*cos(phi)*cos(kfz*zz(k)+phiz)&
                                           +(0.0,1.0)*amp*sin(phi)*sin(kfz*zz(k)+phiz)
            enddo
        endif
        if ((field.eq.'m').or.(field.eq.'b')) then
            do k=1,nlz_par
                smk(j,iloc,k)=smk(j,iloc,k)+amp*cos(phi)*cos(kfz*zz(k)+phiz)&
                                           +(0.0,1.0)*amp*sin(phi)*sin(kfz*zz(k)+phiz)
            enddo
        endif
        !ensure that the field is hermitian when kx=0
        if (i.eq.1) then
            spk(nky-j+2,i,:)=conjg(spk(j,i,:))
            smk(nky-j+2,i,:)=conjg(smk(j,i,:))
        endif
    endif

end subroutine force

subroutine grad(ak,ga)

    use transforms, only: fft,ifft
    use grid, only: kx,ky
    implicit none
    
    real, dimension(:,:,:,:) :: ga
    complex, dimension(:,:,:) :: ak
    complex, allocatable, dimension(:,:,:,:) :: gak
    integer :: i,j,xy

    allocate(gak(nky,nkx_par,nlz_par,2))

    do k=1,nlz_par
        do i=1,nkx_par
            do j=1,nky
                gak(j,i,k,1)=(0.0,1.0)*kx(i)*ak(j,i,k)
                gak(j,i,k,2)=(0.0,1.0)*ky(j)*ak(j,i,k)
            enddo
        enddo
        do xy=1,2
            call ifft(gak(:,:,k,xy),ga(:,:,k,xy))
    enddo

    deallocate(gak)
    deallocate(ak)

end subroutine grad

subroutine crossk(ga1,ga2,bk)
!crosses two real-space 2-vector fields and then goes to k-space
    use transforms, only: fft
    implicit none
    
    real, dimension(:,:,:,:), intent(in) :: ga1,ga2
    complex, dimension(:,:,:),intent(out) :: bk
    real, dimension(:,:,:),allocatable :: b
    integer :: k

    allocate(b(nlx,nly_par,nlz_par))

    b = ga1(:,:,:,1)*ga2(:,:,:,2) - ga1(:,:,:,2)*ga2(:,:,:,1)

    do k=1,nlz_par
        call fft(b(:,:,k),bk(:,:,k))
    enddo

    deallocate(b)

end subroutine crossk

subroutine multkn(ak,k2hak,n)

    use grid, only: kx,ky
    use mp, only: iproc
    implicit none

    real, dimension(:,:,:), intent(in) :: ak
    real, dimension(:,:,:), optional, intent(out) :: k2hak
    integer, optional, intent(in) :: hyper
    integer :: power,k,i,j,iglobal

    if (present(n)) then
        power=n
    else
        power=2
    endif
    
    if (present(k2hak)) then
        do k=1,nlz_par
            do i=1,nkx_par
                iglobal=mod(iproc,npperp)*nkx_par+i
                do j=1,nky
                    !need to take care of 0 mode in case hyper is -ve:
                    !corresponds to removing mean field
                    if ((iglobal.eq.1).and.(j.eq.1) then
                        k2hak(j,i,k)=0.0
                    else
                        k2hak(j,i,k)=(ky(j)**2+kx(i)**2)**power * ak(j,i,k)
                    endif
                enddo
            enddo
        enddo
    else
        do k=1,nlz_par
            do i=1,nkx_par
                iglobal=mod(iproc,npperp)*nkx_par+i
                do j=1,nky
                    if ((iglobal.eq.1).and.(j.eq.1) then
                        ak(j,i,k)=0.0
                    else
                        ak(j,i,k)=(ky(j)**2+kx(i)**2)**power * ak(j,i,k)
                    endif
                enddo
            enddo
        enddo
    endif

end subroutine multkn

subroutine smooth(ak)

    use grid, only: kx,ky
    implicit none

    real, dimension(:,:,:),intent(inout) :: ak
    integer :: i,j,k

    !Hou-Li filtering
    do k=1,nlz_par
        do i=1,nkx_par
            do j=1,nky
                ak(j,i,k)=ak(j,i,k)*exp(-36*( (abs(ky(j)/(nky/2.*lx/ly)))**36 &
                + (abs(kx(i)/((nkx-1)*1.0)))**36))
            enddo
        enddo
    enddo

end subroutine smooth(ak)

real function meanmult(ak,bk)
!   mean of ab using Parseval's theorem
    use mp, only: sum_reduce
    implicit none

    complex, dimension(:,:,:), intent(in) :: ak,bk
    real :: meanmult
    
    meanmult=0.0 
    meanmult=sum(ak(1,:,:)*conjg(bk(1,:,:)))+2.0*sum(a(2,:,:)*conjg(bk(2,:,:)))
    meanmult=meanmult/2.0/nkx/nky/nlz
    call sum_reduce(meanmult,proc0)

end function meanmult

real function rmsk(ak)
!   makes rms of the gradient
    implicit none

    complex, allocatable, dimension(:,:,:) :: dum

    allocate(dum(nky,nkx_par,nlz_par))

    call multkn(ak,dum)
    rmsk=sqrt(meanmult(ak,dum))

end function rmsk


end module subs
