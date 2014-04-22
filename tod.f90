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

integer :: i,j,k,it,ip,irk,isnapfile,ispecfile

real :: time,dt,tlastsnap,tlastspec

real, dimension(:,:,:), allocatable :: zp,zm
complex,dimension(:,:,:), allocatable :: zpk,zmk,spk,smk,dum
real, dimension(:,:,:,2), allocatable :: gzp,gzm,gsp,gsm
logical :: llast=.false.,lout=.false.

character(len=100) :: runname, inputfile
character(len=100) :: filename
character(len=10) :: itstr
character(len=10) :: procstr
character(len=100) :: datadir="data",rundir="./"

!********** Initialization **********

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
call init_mp
call init_grid
call init_transforms

!***** Setting initial fields, if required *****

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
    zp(:,:,:)=0
    zm(:,:,:)=0
endif

!TODO: other initfield options, like "norm"

dt=1e10 !setting timestep to a large number
!alfven speed=1, won't change:
dtva=cfl_frac*dz

if (ladvect) dt=min(dt,dtva)
if (ldiffuse) dt=min(dt,dtv)

!go to fourier space
do k=1,nlz_par
  call fft(zp(:,:,k),zpk(:,:,k))
  call fft(zm(:,:,k),zmk(:,:,k))
enddo

!print out initial values






!start of the timestep
ilastspec=0
ilastsnap=0
tlastspec=0.0
tlastsnap=0.0
it=0
time=0.0
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
    time=time+dt

    !are we going to output to the timeseries?
    lout=(mod(iout,it).eq.0)
    !are we going to stop?
    llast=((it.ge.imax).and.(imax.ge.0)).or.((t.ge.tmax).and.(tmax.ge.0))
    !stopping gracefully - looks for STOP in the run directory
    filename=trim(rundir)//"STOP"
    if (proc0) inquire(file=trim(filename),exist=llast) 
    call broadcast(llast) 

    call smooth(zpk)
    call smooth(zmk)

    zpk0=zpk
    zmk0=zmk

    do irk=rkorder,1,-1
        
        s=dble(irk)
        
        spk=0.0
        smk=0.0
        
        !nl term - reusing names of arrays to save memory
        if (lnonlinear) then
            !gzm and gzp are used twice, so by separating grad from bracket can
            !save some ffts 
            call multkn(zpk,smk)
            call grad(zmk,gzm)
            call grad(smk,gsm)
            call crossk(gzm,gsm,spk) !{zm,k2zp}
            call multkn(zmk,smk)
            call grad(zpk,gzp)
            call grad(smk,gsm)
            call crossk(gzp,gsm,dum) !{zp,k2zm}
            smk=spk+dum
            call multkn(smk,spk,n=-2)
            call crossk(gzp,gzm,dum) !{zp,zm}

            smk=0.5*(spk+dum)
            spk=0.5*(spk-dum)

            call smooth(smk)
            call smooth(spk)

        endif

        if (lout.and.(irk.eq.1)) then
            call multkn(spk,dum)
            work_nlp=meanmult(dum,zp)*dt
            call multkn(smk,dum)
            work_nlm=meanmult(dum,zm)*dt

        !TODO - forcing
        !TODO - work by forcing terms

        !linear (advection) term 
        if (ladvect) then

            call zshift(zpk,+1,dum)
            spk=spk+(znu/dz**2+0.5/dz)*dum
            call zshift(zpk,-1,dum)
            spk=spk+(znu/dz**2-0.5/dz)*dum-2*znu/dz**2*zpk
            
            call zshift(zm,+1,dum)
            smk=smk+(znu/dz**2-0.5/dz)*dum
            call zshift(zm,-1,dum)
            smk=smk+(znu/dz**2+0.5/dz)*dum-2*znu/dz**2*zmk

        endif

        !TODO - work done by this term

        !diffusion term and timestepping
        if (ldiffuse) then

            dum=-nu*dt/s
            call multk2(dum,hyper=hyper_order)
            where (dum.lt.-small)
                dum=exp(dum)
            elsewhere 
                dum=1.0+dum+dum*dum/2.0
            endwhere
            
            zpk=zpk0*dum
            zmk=zmk0*dum

            dum=1-dum
            call divk2(dum,hyper=hyper_order)

            zpk=zpk+spk*dum/nu
            zmk=zmk+smk*dum/nu

        else
            
            !no diffusion
            zpk=zpk0+spk*dt/s
            zmk=zmk0+smk*dt/s

        endif
    
    enddo
   
    !TODO - work done by diffusion

        !TODO - total work

        !TODO - timeseries datafile
        !     - spectra datafile
        !     - snapshots datafile


    !timeseries
    


    !snapshot
    lsnap=((isnap.gt.0).and.(mod(isnap,it).eq.0))
    lsnap=(lsnap.or.((tsnap.gt.0).and.((t-tlastsnap).ge.tsnap)))
    lsnap=(lsnap.or.(llast.and.llastsnap))
    if (lsnap) then
        
        tlastsnap=t
        isnapfile=isnapfile+1 
        write(itstr,"(I0)") isnapfile      
        filename="snap"//trim(itstr)//".dat"
        do k=1,nlz_par    
            zp(:,:,k)=ifft(zpk(:,:,k))
            zm(:,:,k)=ifft(zmk(:,:,k))
        enddo
        call savesnap(filename,zp,zm)
        
    endif
   
    !spectra 
    lspec=((ispec.gt.0).and.(mod(ispec,it).eq.0))
    lspec=(lspec.or.((tspec.gt.0).and.((t-tlastspec).ge.tspec)))
    if (lspec) then
        
        tlastspec=t
        ispecfile=ispecfile+1
        write(itstr,"(I0)") ispecfile
        filename="spec"//trim(itstr)//".dat"
        call spectra(filename,zpk,zmk)
    
    endif
    
enddo timeloop

call finish_mp
end program tod
