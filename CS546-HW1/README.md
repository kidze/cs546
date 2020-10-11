# cs546
cs546 IIT Fall 2020

CS546-HW1 instruction:

In my fusion account (ynguyen), type:

```
cd cs546/CS546-HW1
```

to enter the directory where lies my work for hw1.

Once you are there, you can type:

```
gcc gauss_openmp.c -o gauss_openmp -fopenmp
```

to compile the openmp code, or type:

```
gcc gauss_pthread.c -o gauss_pthread -lpthread
```

to compile the pthread code. Once compile, you can type:

```
./gauss_openmp 10 10 1
```
or

```
./gauss_lpthread 10 10 1 4
```

to run the program. You can add arguments if you want, first and second argument being the Matrix size, third argument being the seed number, and 4th argument being the number of threads desired. The default number of threads I set by default is 4.

Although, in the openmp code I don't use the parameter numThreads so you cannot run the openmp code with 4 arguments.

And there it is! Make sure you play with the arguments so better understand the code. Thank you.
