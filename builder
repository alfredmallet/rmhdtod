ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I.. -c mp_mpi_r8.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I.. -c init.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I.. -c grid.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I.. -c redistribute_mpi.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I.. -c diag.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I.. -c tod.f90
echo "making main program"
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I.. -o tod tod.o mp_mpi_r8.o init.o grid.o redistribute_mpi.o diag.o
