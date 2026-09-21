SingleFile.hpp and SingleFile.cpp contain, respectively, all declarations and the definitions
of the base class and of the exceptions classes defined in Section 4 of the manuscript. 

SingleFileDerived.hpp and SingleFileDerived.cpp contain  declarations and definitions (respectively) of the 3 specific external potentials that are readily implemented in the code. These are, the single file in a constant potential –  SingleFileFlat – , in a linear potential –  SingleFileOnSlope –  and in a harmonic potential –  SingleFileHarmonic. The corresponding analytical single-particle solutions are described in Appendix B in the manuscript.  

IntegerPartitions.hpp and combinations.hpp contain definitions and declarations of these two helper classes detailed in Section 4 and in the Appendix A in the manuscript. 

SingleFileBluePrint.hpp is an example that illustrates how a user should implement his/her own minimal derived class (i.e. a solution for an external potential U(x) that is not readily implemented in the code). When attempting to implement a new potential please follow the instructions in the description of the implementation (Section 4 in the manuscript).     

The attached makefile assumes that the compiler is g++ 7.3.1. The user can adapt it according to his/her favorite compiler.

There are 5 files: figure2.cpp, figure3.cpp, figure4.cpp, figure5a.cpp, figure5b.cpp that generate 5 
different executables (with the same name). Each of them generates (almost; for a limitation of the running-time the naive implementation is omitted when the respective running time becomes long) the plots in the main paper. 

To produce the analogous plots for the remaining two external potential please change in line 23
of the file figure4.cpp: SingleFileHarmonic to SingleFileFlat
or SingleFileOnSlope. 
In line 15 of the same file the variable min=-4 and the variable max=4. This two variables in the case the user is interested in the SingleFileHarmonic define the minimum and the maximum point 
in which the Green's function will be evaluated. If the user is interested in the other two potentials implemented then change the code such that min=0 and max=1, since the solution SingleFileFlat
and SingleFileOnSlope assumes that the reflecting boundary conditions are at this points.

The corresponding running times (CPU: Intel Xeon E3-1270 v2 3.50 GHz 4 cores) are:

1) figure2 runs in around 4 minutes. The naive algorithm is computed only for N<=10 to keep the running time low.
2) figure3 runs in around 3 seconds.
3) figure4 runs in around 13 minutes.
4) figure5a runs in around 6 seconds. The naive algorithm is computed only for N<=10 to keep the running time low.
5) figure5b runs in around 3 minutes.

All the executables use the standard output, except figure4 that generates four .txt files. These are Gnuplot-readable.
