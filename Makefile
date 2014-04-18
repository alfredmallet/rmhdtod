##################################
#        MAKEFILE for TOD        #
##################################


#kraken
FC_kraken=ftn
F90FLAGS_kraken=-target=linux -Mbackslash -r8 -fastsse -O3 -I. -Iutils -I..  -I/opt/fftw/2.1.5.3/include
FLIBS_kraken=-L/opt/fftw/2.1.5.3/lib -ldfftw -ldrfftw

#stampede
FC_stampede = mpif90
F90FLAGS_stampede = -xhost -O2 -r8 -ipo -I${TACC_FFTW2_INC}
FLIBS_stampede = -L${TACC_FFTW2_LIB} -ldrfftw -ldfftw
MKLFLAGS_stampede = -mkl

#set vars
FC = ${FC_${TOD_SYSTEM}}
FLIBS = ${FLIBS_${TOD_SYSTEM}}
F90FLAGS = ${F90FLAGS_${TOD_SYSTEM}}
MKLFLAGS = ${MKLFLAGS_${TOD_SYSTEM}}

OBJS = tod.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o diag.o

.SUFFIXES: .f90

.f90.o:
	$(FC) $(F90FLAGS) -c $<

%.o: %.mod

all: tod

tod:	$(OBJS)
	$(FC) $(F90FLAGS) -o tod $(OBJS) $(FLIBS) $(MKLFLAGS)

debug: F90FLAGS += -g
debug: tod

clean:
	rm -f *.o *.mod

test_make:
	@echo $(USE_FFT)
	@echo $(HOME)

#dependencies
tod.o: init.o mp_mpi_r8.o transforms.o grid.o diag.o
grid.o: init.o fft_work_fftw.o
transforms.o: init.mod fft_work_fftw.o grid.o redistribute_mpi.o
diag.o: mp_mpi_r8.o init.mod grid.o transforms.o
init.o: mp_mpi_r8.o
