program tod
!********************************************************************************
!*                          TOD - solves RMHD equations                         *
!*                                 Alfred Mallet                                *
!*                                                                              *
!*                             based on other codes -                           *
!*            GOSTA(T.Yousef) and VIRIATO(N.Loureiro,A.Kanekar,A.Mallet)        *
!********************************************************************************

use init
use mp
use transforms
use grid
use diag
use subs
implicit none

integer :: i,j,k,it=0,ip,irk

real :: t=0d0,dt,tlastsnap,tlastspec,small,dtva,maxgzx,maxgzy,s,wnlp=0d0,wnlm=0d0
real :: fac,kfz,wfp=0d0,wfm=0d0,wadp=0d0,wadm=0d0,wnup=0d0,wnum=0d0

real, dimension(:,:,:), allocatable :: zp,zm,rdum
complex,dimension(:,:,:), allocatable :: zpk,zmk,spk,smk,dum,dum2,zpk0,zmk0
real, dimension(:,:,:,:), allocatable :: gzp,gzm,gsp,gsm
real, dimension(:,:,:), allocatable :: gzpeq,gzmeq,gg2zpeq,gg2zmeq
logical :: llast=.false.,lout=.false.,lsnap

character(len=100) :: runname, inputfile
character(len=100) :: filename
character(len=10) :: itstr
character(len=10) :: procstr
character(len=100) :: rsloc

integer :: cc,cm
real :: cr

!********** Initialization **********

small=(6*epsilon(1d0)**(1/3.))

call getarg(1,runname)
inputfile=trim(runname)//".in"
call read_parameters(inputfile)

allocate(zp(nlx,nly_par,nlz_par))
allocate(zm(nlx,nly_par,nlz_par))
allocate(gzp(nlx,nly_par,nlz_par,2))
allocate(gzm(nlx,nly_par,nlz_par,2))
allocate(gsp(nlx,nly_par,nlz_par,2))
allocate(gsm(nlx,nly_par,nlz_par,2))

allocate(zpk(nky,nkx_par,nlz_par))
allocate(zmk(nky,nkx_par,nlz_par))
allocate(zpk0(nky,nkx_par,nlz_par))
allocate(zmk0(nky,nkx_par,nlz_par))
allocate(spk(nky,nkx_par,nlz_par))
allocate(smk(nky,nkx_par,nlz_par))
allocate(dum(nky,nkx_par,nlz_par))
allocate(dum2(nky,nkx_par,nlz_par))

allocate(rdum(nky,nkx_par,nlz_par))

if (lequil) then
    allocate(gzpeq(nlx,nly_par,2))
    allocate(gzmeq(nlx,nly_par,2))
    allocate(gg2zpeq(nlx,nly_par,2))
    allocate(gg2zmeq(nlx,nly_par,2))
endif

call init_mp
call init_grid
call init_transforms

call system_clock(cc,cr,cm)
cr=uniran(init=cc)

!*****Loading equilibrium field, if required. *****
!Also getting pieces that are needed in the eq-pert interaction term.
!TODO: check dummies are reset properly!
if (lequil) then
    if (equilmode="func") then
        call makeeq(equilfunc,zp(:,:,1),zm(:,:,1))
    else
        call loadeq(equilfile,equilpath,zp(:,:,1),zm(:,:,1))
    call fft(zp(:,:,1),zpk(:,:,1))
    call grad2d(zpk(:,:,1),gzpeq)
    call fft(zm(:,:,1),zmk(:,:,1))
    call grad2d(zmk(:,:,1),gzmeq)
    call multkn2d(zpk(:,:,1),dum(:,:,1))
    call grad2d(dum(:,:,1),gg2zpeq)
    call multkn2d(zmk(:,:,1),dum(:,:,1))
    call grad2d(dum(:,:,1),gg2zmeq)
endif


!***** Setting initial fields, if required *****
if (restart.eq.0) then
    if (initfield.eq."wave") then
        if (proc0) then
            write(*,*) "Sinusoidal initial field with given wavenumbers and amplitudes, kip=",kipx,kipy,kipz,"kim=",kimx,kimy,kimz,"ampzp=",ampzp,"ampzm=",ampzm
        endif
        do k=1,nlz_par
            do j=1,nly_par
                do i=1,nlx
                    zp(i,j,k)=ampzp*sin(2.*pi*(kipx*xx(i)/lx+kipy*yy(j)/ly+kipz*zz(k)/lz))
                    zm(i,j,k)=ampzm*sin(2.*pi*(kimx*xx(i)/lx+kimy*yy(j)/ly+kimz*zz(k)/lz))
                enddo
            enddo
        enddo
    else
        zp(:,:,:)=0d0
        zm(:,:,:)=0d0
    endif
    isnapfile=0
    tlastsnap=0.0
    it=0
    t=0d0
else
    call loadsnap(rspath,rsfile,zp,zm)
    isnapfile=irsfile
    tlastsnap=rstime
    it=itrs
    t=rstime
endif 
call barrier
filename="start.dat"
call savesnap(filename,zp,zm,t)

!TODO: other initfield options, like "norm"

dt=1d10 !setting timestep to a large number
!alfven speed=1, won't change:
dtva=cfl_frac*dz

if (ladvect) dt=min(dt,dtva)

!go to fourier space
do k=1,nlz_par
  call fft(zp(:,:,k),zpk(:,:,k))
  call fft(zm(:,:,k),zmk(:,:,k))
enddo

!print out initial values

call outputts(tsfile)

!start of the timestep
!isnapfile=0
!tlastsnap=0.0
!it=0
!t=0.0

timeloop: do

    !max value of perp flows
    call grad(zpk,gzp)
    call grad(zmk,gzm)
    maxgzx=max(maxval(abs(gzp(:,:,:,2))),maxval(abs(gzm(:,:,:,2))))
    maxgzy=max(maxval(abs(gzp(:,:,:,1))),maxval(abs(gzm(:,:,:,1))))
    call max_allreduce(maxgzx)
    call max_allreduce(maxgzy)

    if (lnonlinear) dt=min(dt,cfl_frac*dx/maxgzx,cfl_frac*dy/maxgzy)

    !update time and iteration no.
    it=it+1
    t=t+dt

    !are we going to output to the timeseries?
    lout=(mod(it,iout).eq.0)
    !stopping gracefully - looks for STOP in the run directory
    filename=trim(rundir)//"STOP"
    if (proc0) inquire(file=trim(filename),exist=llast) 
    call broadcast(llast) 
    !are we going to stop?
    llast=((llast).or.((it.ge.imax).and.(imax.ge.0)).or.((t.ge.tmax).and.(tmax.ge.0)))

    call smooth(zpk)
    call smooth(zmk)

    zpk0=zpk
    zmk0=zmk

    do irk=rkorder,1,-1
        
        s=dble(irk)
        
        spk=0d0
        smk=0d0
        
        !nl term - reusing names of arrays to save memory
        if (lnonlinear) then
            !gzm and gzp are used twice, so by separating grad from bracket can
            !save some ffts
            !saving one fft by changing to cross instead of crossk
            call multkn(zpk,smk)
            call grad(zmk,gzm)
            call grad(smk,gsm)
            call cross(gzm,gsm,gsp) !{zm,k2zp}
            !call crossk(gzm,gsm,spk) !{zm,k2zp}
            call multkn(zmk,smk)
            call grad(zpk,gzp)
            call grad(smk,gsm)
            call cross(gzp,gsm,gzm) !{zp,k2zm}
            !call crossk(gzp,gsm,dum) !{zp,k2zm}
            !smk=spk+dum
            gsp=gsp+gzm
            do k=1,nlz_par
                call fft(gsp(:,:,k),smk(:,:,k))
            enddo
            call multkn(smk,spk,n=-2)
            call crossk(gzp,gzm,dum) !{zp,zm}

            smk=-0.5d0*(spk+dum)
            spk=-0.5d0*(spk-dum)

            call smooth(smk)
            call smooth(spk)

        endif
        
        ! work done by nonlinear term
        if (lout.and.(irk.eq.1)) then
            call multkn(spk,dum)
            wnlp=meanmult(dum,zpk)*dt
            call multkn(smk,dum)
            wnlm=meanmult(dum,zmk)*dt
        endif

        ! forcing
        if (lforce) then 
            if (epsp.ge.0) then
                fac=lz/2./pi
                kfz=nint(kfz1*fac-0.5*uniran()*(fac*kfz2-fac*kfz1+1.0))/fac
                kfz=kfz*nint(sign(1.0,(2.0*uniran()-1.0)))
                if (epsm.ge.0) then !elsasser forcing
                    if (epsp.gt.0) call force(kfz,epsp,dt/s,spk,smk,'p')
                    kfz=nint(kfz1*fac-0.5*uniran()*(fac*kfz2-fac*kfz1+1.0))/fac
                    kfz=kfz*nint(sign(1.0,(2.0*uniran()-1.0)))
                    if (epsm.gt.0) call force(kfz,epsm,dt/s,spk,smk,'m')
                else
                    if (epsp.gt.0) call force(kfz,epsp,dt/s,spk,smk,'b')
                endif
            endif
        endif

        ! work done by forcing terms
        if (lout.and.(irk.eq.1)) then
            call multkn(spk,dum)
            wfp=meanmult(dum,zpk)*dt-wnlp
            call multkn(smk,dum)
            wfm=meanmult(dum,zmk)*dt-wnlm 
        endif

        !linear (advection) term 
        if (ladvect) then
            call zshift(zpk,+1,dum)
            spk=spk+(znu/dz**2+0.5/dz)*dum
            call zshift(zpk,-1,dum)
            spk=spk+(znu/dz**2-0.5/dz)*dum-2*znu/dz**2*zpk
            
            call zshift(zmk,+1,dum)
            smk=smk+(znu/dz**2-0.5/dz)*dum
            call zshift(zmk,-1,dum)
            smk=smk+(znu/dz**2+0.5/dz)*dum-2*znu/dz**2*zmk

        endif

        !work done by advection term
        if (lout.and.(irk.eq.1)) then
            call multkn(spk,dum)
            wadp=meanmult(dum,zpk)*dt-wnlp-wfp
            wadm=meanmult(dum,zmk)*dt-wnlm-wfm
        endif

        !diffusion term and timestepping
        if (ldiffuse) then
            dum=-nu*dt/(1d0*s)
            call multkn(dum,n=2*hyper_order)
            where (abs(dum).gt.small)
                dum2=exp(-abs(dum))
                dum=1d0-dum2
            elsewhere 
                dum2=abs(1d0+dum+dum*dum/2d0)
                dum=-dum-dum*dum/2d0
            endwhere
            zpk=zpk0*dum2
            zmk=zmk0*dum2
            call multkn(dum,n=-2*hyper_order)
            zpk=zpk+spk*dum/nu
            zmk=zmk+smk*dum/nu
            if (mod(iproc,npperp).eq.0) then 
                zpk(1,1,:)=0d0
                zmk(1,1,:)=0d0
            endif
        
        else
            
            !no diffusion
            zpk=zpk0+spk*dt/s
            zmk=zmk0+smk*dt/s

        endif
    
    enddo
   
    !work done by diffusion
    if (lout) then
        call multkn(zpk0,dum,n=hyper_order+1)
        wnup=meanmult(dum,dum)
        call multkn(zmk0,dum,n=hyper_order+1)
        wnum=meanmult(dum,dum)
    endif
        
    !timeseries 
    if (lout) call outputts(tsfile)

    !snapshot
    lsnap=((isnap.gt.0).and.(modulo(it,isnap).eq.0))
    lsnap=(lsnap.or.((tsnap.gt.0).and.((t-tlastsnap).ge.tsnap)))
    lsnap=(lsnap.or.(llast.and.llastsnap))
    if (lsnap) then
        
        tlastsnap=t
        isnapfile=isnapfile+1 
        write(itstr,"(I4.4)") isnapfile      
        filename="snap"//trim(itstr)//".dat"
        do k=1,nlz_par    
            call ifft(zpk(:,:,k),zp(:,:,k))
            call ifft(zmk(:,:,k),zm(:,:,k))
        enddo
        call savesnap(filename,zp,zm,t)
        
    endif
   
    !spectra - TODO: write subroutine
    !lspec=((ispec.gt.0).and.(mod(ispec,it).eq.0))
    !lspec=(lspec.or.((tspec.gt.0).and.((t-tlastspec).ge.tspec)))
    !if (lspec) then
    !    
    !    tlastspec=t
    !    ispecfile=ispecfile+1
    !    write(itstr,"(I0)") ispecfile
    !    filename="spec"//trim(itstr)//".dat"
    !    call spectra(filename,zpk,zmk,t)
    !
    !endif
    if (llast) exit timeloop

enddo timeloop

call finish_mp

contains

subroutine outputts(filename)
    
    implicit none
    
    character (len=100),optional :: filename
    logical, save :: lfirst=.true.
    complex, dimension(:,:,:),allocatable :: dum,dum2,phik,psik
    real :: ep,em,eu,eb,er,h,rmszp,rmszm
    character(len=1000) :: astring

    allocate(dum(nky,nkx_par,nlz_par))
    allocate(dum2(nky,nkx_par,nlz_par))
    allocate(phik(nky,nkx_par,nlz_par))
    allocate(psik(nly,nkx_par,nlz_par))

    !energies and amplitudes
    call multkn(zpk,dum2)
    ep=meanmult(zpk,dum2) ! energy in z+ fluctuations
    call multkn(zmk,dum)    
    em=meanmult(zmk,dum) ! energy in z- fluctuations
    phik=0.5*(zpk+zmk)
    call multkn(phik,dum)
    eu=meanmult(phik,dum) ! energy in velocity fluctuations
    psik=0.5*(zpk-zmk)
    call multkn(psik,dum)
    eb=meanmult(psik,dum) ! energy in magnetic fluctuations
    rmszp=sqrt(meanmult(zpk,zpk)) ! rms z+
    rmszm=sqrt(meanmult(zmk,zmk)) ! rms z-
    er=meanmult(zmk,dum2) ! <z+.z->, residual energy
    h=meanmult(phik,dum) ! <u.b>, cross-helicity

    if (proc0.and.lfirst) then
        open(35,file=trim(datadir)//"/"//trim(filename),recl=10000,position="append")
        write(astring,'(A10,18A25)') &
                       "    it    ","            t            ",&
        "           Ep            ","           Em            ",&
        "           Eu            ","           Eb            ",&
        "           Er            ","           H             ",&
        "           zp            ","           zm            ",&
        "          Wnup           ","          Wnum           ",&
        "          Wnlp           ","          Wnlm           ",&
        "          Wadp           ","          Wadm           ",& 
        "           Wfp           ","           Wfm           ",&
        "           dt            "
        write(35,*) trim(astring)
        write(*,*) trim(astring)
        lfirst=.false.
        close(35)
    endif

    if (proc0) then 
        write(astring,'(I10,18G25.8)') it,t,ep,em,eu,eb,er,h,rmszp,rmszm,wnup,wnum,wnlp,wnlm,wadp,wadm,wfp,wfm,dt
        !write(astring,'(I10,18G25.8)') it,t,Ep,Em,Eu,Eb,Er,H,zp,zm,wnup,wnum,wnlp,wnlm,wadp,wadm,wfp,wfm,dt
        if (present(filename)) then
            open(35,file=trim(datadir)//"/"//trim(filename),recl=10000,position="append")
            write(35,*) trim(astring)
            close(35)
        endif
        write(*,*) trim(astring)
    endif

    deallocate(dum)
    deallocate(dum2)
    deallocate(phik)
    deallocate(psik)

end subroutine outputts

end program tod
