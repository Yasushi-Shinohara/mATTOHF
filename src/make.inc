#FC options
FC = mpif90 -fopenmp #With GNU



#LN options
LN = -L/usr/lib64/ -llapack -lblas  -I/usr/include/ -lfftw3 #For OS-default libraries in my local environment

