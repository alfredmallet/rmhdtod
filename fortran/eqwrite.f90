program eqwrite

!******* Writes an equilibrium into a file **********

use init
use mp
use transforms
use grid
use diag
use subs
implicit none

call getarg(1,runname)
inputfile=trim(runname)//".in"
call read_parameters(inputfile

!allocate stuff

allocate(zp(nlx,nly_par))
allocate(zm(nlx,nly_par))

call init_mp
call init_grid

!
