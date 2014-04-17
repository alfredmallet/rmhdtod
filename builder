ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -c mp_mpi_r8.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -c init.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -c fft_work_fftw.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -c grid.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -c redistribute_mpi.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -c transforms.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -c diag.f90
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -c tod.f90
echo "Compiling main program"
ftn -target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include   -o tod tod.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o diag.o  -L/opt/fftw/2.1.5.3/lib -ldfftw -ldrfftw
