# aux2psd
This code reads the WRF SBM auxiliary files and convert them to PSD related variables.

To compile the fortran conversion files, one needs to:
* gfortran -c -fPIC readSBM.f90
* f2py -c -m aux2py aux2py.f90 readSBM.o

These steps create a python module named aux2py