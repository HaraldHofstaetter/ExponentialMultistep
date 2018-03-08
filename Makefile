INC_FFTW3 = /usr/include
LIB_FFTW3 = /usr/lib

#INC_FFTW3 = $(HOME)/include
#LIB_FFTW3 = $(HOME)/lib

EXECUTABLES = test_gaussian test_soliton

all: $(EXECUTABLES)

phi_functions.o phi_functions.mod: phi_functions.f90
	gfortran -g -c -o phi_functions.o phi_functions.f90

wavefunctions_fourier1d.o wavefunctions_fourier1d.mod:  phi_functions.mod wavefunctions_fourier1d.f90
	gfortran -g -c -I$(INC_FFTW3) -o wavefunctions_fourier1d.o wavefunctions_fourier1d.f90

test_gaussian.o: wavefunctions_fourier1d.mod test_gaussian.f90
	gfortran -g -c -o test_gaussian.o test_gaussian.f90

test_gaussian: phi_functions.o wavefunctions_fourier1d.o test_gaussian.o
	gfortran -g -o test_gaussian test_gaussian.o \
	    wavefunctions_fourier1d.o  phi_functions.o \
	    -L$(LIB_FFTW3) -lfftw3

exponential_multistep.o exponential_multistep.mod: wavefunctions_fourier1d.mod exponential_multistep.f90
	gfortran -g -c -o exponential_multistep.o exponential_multistep.f90

test_soliton.o: exponential_multistep.mod test_soliton.f90
	gfortran -g -c -o test_soliton.o test_soliton.f90

test_soliton: phi_functions.o wavefunctions_fourier1d.o exponential_multistep.o test_soliton.o
	gfortran -g -o test_soliton test_soliton.o \
	    exponential_multistep.o wavefunctions_fourier1d.o phi_functions.o \
	    -L$(LIB_FFTW3) -lfftw3

clean:
	rm *.o *.mod $(EXECUTABLES)
