HOW TO COMPILE AND RUN CODE

locate to the cs546 file in ynguyen
locate CS546-HW3 folder

inside there should be a file name
gauss_mpi.c

type this to compile

mpicc -o gauss_mpi gauss_mpi.c

then type this to run the file

mpiexec -np <number of processors> ./gauss_mpi <matrix size>


for example run this:

mpiexec -np 4 ./gauss_mpi 3000

4 processors and N=3000 sized matrix

Observe the elapsed time. Thank you!
