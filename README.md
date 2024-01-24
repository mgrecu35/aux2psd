# aux2psd
read SBM auxiliary files

To compile the fortran conversion files, one needs to:
* gfortran -c -fPIC readSBM.f90
* f2py -c -m aux2py aux2py.f90 readSBM.o

These steps create a python module named aux2py