# star_discr_subset
Branch and bound algorithm to solve the star discrepancy subset selection problem

The folder "Simple Search" contains a recursive algorithm that evaluates all solutions.
The folder "Branch and Bound" contains two branch and bound algorithms to solve this problem in 2D (main_sort.cpp) and for any number of dimensions (main_bound.cpp)

To compile the program use the command:
g++ -O3 [file].cpp *.c -o main

And then run the program with the command:
./main [k] < [input file]

The files bzdiscr.c, bzdiscr.h, dem_discr.c and dem_discr.h were coded by Magnus Wahlström and where then adapted to be compilable by g++
