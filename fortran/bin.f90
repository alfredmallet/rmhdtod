program binner

!Bins sf produced by sf3.f90. Designed to be run on 32cpu

use init
use mp
use transforms
use grid
use diag
use subs
implicit none

character (len=100) :: runname, inputfile, sfdir,wfile,procstr

integer :: isep, nsep=32,isamp,nsamp=500000,iang,iali,nang=9,im,nm=10,wth

real :: thetaBloc,thetapm,thetaub
real, dimension(4) :: amps, thetas
real, dimension(:),allocatable :: thetabins
real, dimension(:,:,:,:,:,:),allocatable :: sf
real, dimension(:,:,:,:,:),allocatable :: lsf
integer, dimension(:,:,:,:),allocatable :: ninbin
integer, dimension(:), allocatable :: nperp
real, dimension(:,:), allocatable :: thetapmn,thetaubn,botpm,botub

allocate(sf(4,nm,nang,4,nang,nsep))
allocate(ninbin(nang,4,nang,nsep))
allocate(thetabins(nang))
allocate(thetapmn(nm,nsep))
allocate(thetaubn(nm,nsep))
allocate(nperp(nm))
allocate(botpm(nm,nsep))
allocate(botub(nm,nsep))
allocate(lsf(4,nang,4,nang,nsep))

call getarg(1,runname)
inputfile=trim(runname)//".in"
call read_parameters(inputfile)
call getarg(2,sfdir)

call init_mp

write(procstr,"(I0)") iproc

do isep=1,nsep
    write(wfile,"(I0)") isep
    open(10,file=trim(sfdir)//"/"//trim(procstr)//"/base"//trim(wfile)//".dat")

    do isamp=1,nsamp
        
        read(10,"(11G25.8)") amps(1),amps(2),amps(3),amps(4),thetaBloc,thetas(1),thetas(2),thetas(3),thetas(4),thetapm,thetaub
        if (floor(asin(sin(thetaBloc))*2*nang/pi+1).eq.nang) then
            nperp(isep)=nperp(isep)+1
            do im=1,nm
                thetapmn(im,isep)=thetapmn(im,isep)+(amps(1)*amps(2)*abs(sin(thetapm)))**(0.5*im)
                thetaubn(im,isep)=thetaubn(im,isep)+(amps(3)*amps(4)*abs(sin(thetaub)))**(0.5*im)
                botpm(im,isep)=botpm(im,isep)+(amps(1)*amps(2))**(0.5*im)
                botub(im,isep)=botub(im,isep)+(amps(1)*amps(2))**(0.5*im)
            enddo
        endif
        if (amps(1)>100.) write(*,*) isep, iproc, isamp
        do iang=1,nang
            if (floor(asin(sin(thetaBloc))*2*nang/pi+1).eq.iang) then
                do wth=1,4
                    do iali=1,nang
                        if (floor(asin(sin(thetas(wth)))*2*nang/pi+1).eq.iali) then
                            ninbin(iali,wth,iang,isep)=ninbin(iali,wth,iang,isep)+1
                            lsf(:,iali,wth,iang,isep)=lsf(:,iali,wth,iang,isep)+log(amps(:))
                            do im=1,nm
                                sf(:,im,iali,wth,iang,isep)=sf(:,im,iali,wth,iang,isep)+amps(:)**(im*0.5) 
                            enddo 
                        endif
                    enddo
                enddo
            endif
         enddo
    enddo
    
enddo
thetapmn=thetapmn/botpm
thetaubn=thetaubn/botub
call barrier
do isep=1,nsep
do im=1,nm
call sum_reduce(thetapmn(im,isep),0)
call sum_reduce(thetaubn(im,isep),0)
enddo
do iang=1,nang
do wth=1,4
do iali=1,nang
call sum_reduce(ninbin(iali,wth,iang,isep),0)
call sum_reduce(lsf(:,iali,wth,iang,isep),0)
do im=1,nm
call sum_reduce(sf(:,im,iali,wth,iang,isep),0)
enddo
enddo
enddo
enddo
enddo
call barrier
write(*,*) "gothere", iproc
if (iproc.eq.0) then
    do wth=1,4
        lsf(wth,:,:,:,:)=lsf(wth,:,:,:,:)/(ninbin(:,:,:,:)*1.0)
    do im=1,nm
        sf(wth,im,:,:,:,:)=sf(wth,im,:,:,:,:)/(ninbin(:,:,:,:)*1.0)
    enddo
    enddo
    open(15,file=trim(sfdir)//"/3dsf.dat")
    open(17,file=trim(sfdir)//"/lsf.dat")
    do isep=1,nsep
    do iang=1,nang
    do wth=1,4
    do iali=1,nang
        write(17,"(1G25.8)") lsf(1,iali,wth,iang,isep),lsf(2,iali,wth,iang,isep),lsf(3,iali,wth,iang,isep),lsf(4,iali,wth,iang,isep)
    do im=1,nm
        write(15,"(10G25.8)") isep, iang, wth, iali, im, ninbin(iali,wth,iang,isep), sf(1,im,iali,wth,iang,isep),sf(2,im,iali,wth,iang,isep),sf(3,im,iali,wth,iang,isep),sf(4,im,iali,wth,iang,isep)
    enddo
    enddo
    enddo
    enddo
    enddo
    close(15)
    open(16,file=trim(sfdir)//"/angsf.dat")
    do isep=1,nsep
    do im=1,nm
        write(16,"(4G25.8)") isep,thetapmn(im,isep),thetaubn(im,isep),nperp(isep)
    enddo
    enddo
    close(16)
    close(17)
    write(*,*) 'im done'
endif
call finish_mp
end program binner
