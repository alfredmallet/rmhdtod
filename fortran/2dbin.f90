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

integer :: isep, nsep=32,isamp,nsamp=256000,iali,nang=9,im,nm=20,wth

real :: thetapm,thetaub
real, dimension(4) :: amps, thetas
real, dimension(:),allocatable :: thetabins
real, dimension(:,:,:,:,:),allocatable :: sf
integer, dimension(:,:,:),allocatable :: ninbin

allocate(sf(4,nm,nang,4,nsep))
allocate(ninbin(nang,4,nsep))
allocate(thetabins(nang))

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
        
        read(10,"(10G25.8)") amps(1),amps(2),amps(3),amps(4),thetas(1),thetas(2),thetas(3),thetas(4),thetapm,thetaub
        if (amps(1)>100.) write(*,*) isep, iproc, isamp
                do wth=1,4
                    do iali=1,nang
                        if (floor(asin(sin(thetas(wth)))*2*nang/pi+1).eq.iali) then
                            ninbin(iali,wth,isep)=ninbin(iali,wth,isep)+1
                            do im=1,nm
                                sf(:,im,iali,wth,isep)=sf(:,im,iali,wth,isep)+amps(:)**(im*0.5) 
                            enddo 
                        endif
                    enddo
                enddo
    enddo
    

enddo
call barrier
do isep=1,nsep
do wth=1,4
do iali=1,nang
call sum_reduce(ninbin(iali,wth,isep),0)
do im=1,nm
call sum_reduce(sf(:,im,iali,wth,isep),0)
enddo
enddo
enddo
enddo
call barrier
if (iproc.eq.0) then
    do wth=1,4
    do im=1,nm
        sf(wth,im,:,:,:)=sf(wth,im,:,:,:)/(ninbin(:,:,:)*1.0)
    enddo
    enddo
    open(15,file=trim(sfdir)//"/3dsf.dat")
    do isep=1,nsep
    do wth=1,4
    do iali=1,nang
    do im=1,nm
        write(15,"(10G25.8)") isep, wth, iali, im, ninbin(iali,wth,isep), sf(1,im,iali,wth,isep),sf(2,im,iali,wth,isep),sf(3,im,iali,wth,isep),sf(4,im,iali,wth,isep)
    enddo
    enddo
    enddo
    enddo
endif

call finish_mp
end program binner
