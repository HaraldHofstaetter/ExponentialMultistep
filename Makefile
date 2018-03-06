INC_FFTW3 = /usr/include
LIB_FFTW3 = /usr/lib

wavefunctions_fourier1d.o: wavefunctions_fourier1d.f90
	gfortran -g -c -I$(INC_FFTW3) -o wavefunctions_fourier1d.o wavefunctions_fourier1d.f90

test_gaussian.o: test_gaussian.f90
	gfortran -g -c -o test_gaussian.o test_gaussian.f90

test_gaussian: test_gaussian.o wavefunctions_fourier1d.o
	gfortran -g -o test_gaussian test_gaussian.o wavefunctions_fourier1d.o -L$(LIB_FFTW3) -lfftw3

exponential_multistep.o: exponential_multistep.f90
	gfortran -g -c -o exponential_multistep.o exponential_multistep.f90
