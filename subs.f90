module subs

use init

implicit none

contains

subroutine force(kfz,eps,dt,spk,smk,field)

    use init
    use mp, only: proc0,iproc
    use grid, only: zz
    implicit none

    real, intent(in) :: dt,eps,kfz
    character, intent(in) :: field
    integer, dimension(:,:), allocatable :: kfs_big
    integer, dimension(:,:), allocatable, save :: kfs
    complex, dimension(:,:,:) :: spk,smk
    real :: kfp, kav=0
    real :: amp, phi, phiz
    integer :: i,j,ik,k,nk=0,sgn,proc_force,iloc
    logical, save :: lfirst=.true.

    ! calculate possible elements
    if (lfirst) then
        allocate(kfs_big(nkx*nky,2))
        do i=-ceiling(kfp2),ceiling(kfp2)
            do j=-ceiling(kfp2),ceiling(kfp2)
                kfp=sqrt((1.0*i)**2+(1.0*j)**2)
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
    amp=sqrt(-eps/dt/kfp**2.0*log(uniran()))*sqrt(1.0*nlx*nly) !N.B. may depend on fft convention
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
        if ((field.eq.'b').or.(field.eq.'p')) then
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
    integer :: i,j,k,xy

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
    enddo

    deallocate(gak)

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

subroutine multkn(ak,knhak,n)

    use grid, only: kx,ky
    use mp, only: iproc
    implicit none

    complex, dimension(:,:,:), intent(inout) :: ak
    complex, dimension(:,:,:), optional, intent(out) :: knhak
    integer, optional, intent(in) :: n
    integer :: power,k,i,j,iglobal

    if (present(n)) then
        power=n
    else
        power=2
    endif
    
    if (present(knhak)) then
        do k=1,nlz_par
            do i=1,nkx_par
                iglobal=mod(iproc,npperp)*nkx_par+i
                do j=1,nky
                    !need to take care of 0 mode in case n is -ve:
                    !corresponds to removing mean field
                    if ((iglobal.eq.1).and.(j.eq.1)) then
                        knhak(j,i,k)=0.0
                    else
                        knhak(j,i,k)=(ky(j)**2+kx(i)**2)**power * ak(j,i,k)
                    endif
                enddo
            enddo
        enddo
    else
        do k=1,nlz_par
            do i=1,nkx_par
                iglobal=mod(iproc,npperp)*nkx_par+i
                do j=1,nky
                    if ((iglobal.eq.1).and.(j.eq.1)) then
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

    complex, dimension(:,:,:),intent(inout) :: ak
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

end subroutine smooth

real function meanmult(ak,bk)
!   mean of ab using Parseval's theorem
    use mp, only: sum_reduce
    implicit none

    complex, dimension(:,:,:), intent(in) :: ak,bk
    
    meanmult=0.0 
    meanmult=sum(ak(1,:,:)*conjg(bk(1,:,:)))+2.0*sum(ak(2,:,:)*conjg(bk(2,:,:)))
    meanmult=meanmult/2.0/nkx/nky/nlz
    call sum_reduce(meanmult,0)

end function meanmult

real function rmsk(ak)
!   makes rms of the gradient
    implicit none

    complex, allocatable, dimension(:,:,:) :: dum
    complex, dimension(:,:,:) :: ak

    allocate(dum(nky,nkx_par,nlz_par))

    call multkn(ak,dum)
    rmsk=sqrt(meanmult(ak,dum))

end function rmsk

function uniran(init)
!
! 26-sep-02/wolf: Adapted from `Numerical Recipes for F90' ran() routine
! 23-mai-06/tay: Ripped from Pencil-code
! 24/4/14: AM ripped from gosta
!
! "Minimal" random number generator of Park and Miller combined
! with a Marsaglia shift sequence. Returns a uniform random deviate
! between 0.0 and 1.0 (exclusive of the endpoint values).
! Call with (INIT=ival) to initialize.
! The period of this generator is supposed to be about 3.1× 10^18.
!
implicit none
!
real (kind=8) :: uniran
integer, parameter :: mseed=256
! integer, dimension(mseed), save :: seed=0
integer, dimension(mseed), save :: rstate=0
real (kind=8), parameter :: impossible=3.9085e37
real (kind=8), save :: am=impossible  ! will be constant on a given platform
integer, optional, intent(in) :: init
integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
integer :: k,init_ts=1812   ! default value
logical, save :: first_call=.true.

!ajw This doesn't appear to always get set!
if (first_call) then
  am=nearest(1.0,-1.0)/im
  first_call=.false.
endif
if (present(init) .or. rstate(1)==0 .or. rstate(2)<=0) then
  !
  ! initialize
  !
  if (present(init)) init_ts = init
  am=nearest(1.0,-1.0)/im
  rstate(1)=ieor(777755555,abs(init_ts))
  rstate(2)=ior(ieor(888889999,abs(init_ts)),1)
endif
!
! Marsaglia shift sequence with period 2^32-1
!
rstate(1)=ieor(rstate(1),ishft(rstate(1),13))
rstate(1)=ieor(rstate(1),ishft(rstate(1),-17))
rstate(1)=ieor(rstate(1),ishft(rstate(1),5))
!
! Park-Miller sequence by Schrage's method, period 2^31-2
!
k=rstate(2)/iq
rstate(2)=ia*(rstate(2)-k*iq)-ir*k
if (rstate(2) < 0) rstate(2)=rstate(2)+im
!
! combine the two generators with masking to ensure nonzero value
!
uniran=am*ior(iand(im,ieor(rstate(1),rstate(2))),1)
endfunction uniran

end module subs
