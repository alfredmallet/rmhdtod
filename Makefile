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

OBJS = tod.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o subs.o diag.o 

VECOBJS = vecwrite.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o subs.o diag.o

SFOBJS = sf.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o subs.o diag.o

SF2OBJS = sf2d.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o subs.o diag.o

SF3OBJS = sf3d.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o subs.o diag.o

BINOBJS = bin.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o subs.o diag.o

RCBOBJS = rcb.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o subs.o diag.o

2DBINOBJS = 2dbin.o fft_work_fftw.o init.o transforms.o grid.o mp_mpi_r8.o redistribute_mpi.o subs.o diag.o

.SUFFIXES: .f90

.f90.o:
	$(FC) $(F90FLAGS) -c $<

%.o : %.mod

all: tod

tod:	$(OBJS)
	$(FC) $(F90FLAGS) -o tod $(OBJS) $(FLIBS) $(MKLFLAGS)

debug: F90FLAGS += -g
debug: tod

vec: vecwrite

vecwrite: $(VECOBJS)
	  $(FC) $(F90FLAGS) -o vecwrite $(VECOBJS) $(FLIBS) $(MKLFLAGS)

sf2d:	$(SF2OBJS)
	$(FC) $(F90FLAGS) -o sf2d $(SF2OBJS) $(FLIBS) $(MKLFLAGS)


sf3d:	$(SF3OBJS)
	$(FC) $(F90FLAGS) -o sf3d $(SF3OBJS) $(FLIBS) $(MKLFLAGS)

bin:	$(BINOBJS)
	$(FC) $(F90FLAGS) -o bin $(BINOBJS) $(FLIBS) $(MKLFLAGS)

2dbin: $(2DBINOBJS)
	$(FC) $(F90FLAGS) -o 2dbin $(2DBINOBJS) $(FLIBS) $(MKLFLAGS)

rcb:	$(RCBOBJS)
	$(FC) $(F90FLAGS) -o rcb $(RCBOBJS) $(FLIBS) $(MKLFLAGS)

clean:
	rm -f *.o *.mod

test_make:
	@echo $(USE_FFT)
	@echo $(HOME)
	@echo $(FC)
	@echo $(F90FLAGS)
	@echo $(FC_stampede)
	@echo ${F90FLAGS_${TOD_SYSTEM}}
	@echo ${TOD_SYSTEM}
#dependencies
tod.o: init.o mp_mpi_r8.o transforms.o grid.o diag.o subs.o
vecwrite.o: init.o mp_mpi_r8.o transforms.o grid.o diag.o subs.o
grid.o: init.o fft_work_fftw.o
transforms.o: init.mod fft_work_fftw.o grid.o redistribute_mpi.o
diag.o: mp_mpi_r8.o init.mod grid.o transforms.o
subs.o: mp_mpi_r8.o init.mod grid.o transforms.o
init.o: 
mp_mpi_r8.o: init.mod
