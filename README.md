Documentation                                               {#mainpage}
=============

Version
-------
cmc++ 1.0.0 - Sep 2016.

What is new
-----------
Start version ready.

Download
--------
If you are in the user list of the project's [Bitbucket repository][bb], do a

	hg clone https://amedhi@bitbucket.org/amedhi/classical-mc 

You can also click [here][dl] to download. Otherwise send me a personal request 
at 'amedhi@iisertvm.ac.in'.

[bb]: https://bitbucket.org/amedhi/classical-mc
[dl]: https://bitbucket.org/amedhi/classical-mc/downloads

Introduction
------------

cmc++ is a program written in C++ for doing Monte Carlo simulation for 
'classical' model Hamiltonians like the Ising model. 
The core codes are model and lattice geometry independent. This means 
it should work for all 'classical' models and all lattice geometries once 
these definitions are incorporated into the project in an appropriate manner. 
The program can be build in serial mode or in MPI mode by simply a choice of
compiler directive. In parallel mode, the program creates multiple 
independent Markov chains in each processor and the final results 
averaged over the results from these independent chains.  The code is designed 
to make it easy to add new features - such as, adding a new model or lattice, 
adding a new physical observable to calculate or a new input parameter.

Installation
------------

There are two dependencies - Boost and Eigen C++. You will also need 
latest C++ compilers as the code contains many C++11 features. It works with
latest g++ and clang++. It should work with others also.

To build, edit the included Makefile manually to specify your system 
variables. To build the MPI version, use the preprocessor directive 
'HAVE_MPI' with an appropriate compiler as suggested in the Makefile. 
Then do

	make 

To install the executable, do

	make install 


[mpi]: http://www.open-mpi.org/

Usage
-----
Once the executable is built and installed, run it as 

	program [OPTIONS] [FILE]  

Here 'program' is the name of the executable.  The program accepts a few 
optional arguments OPTIONS. Type 

	program -help

to see the list of available options. The FILE argument, if present, 
is taken as the input filename, else a default filename 'input.parm' is 
looked for in the current directory.  The format of the input file along 
with the descriptions of the input parameters is given in the 
\ref inputdoc. An example command to run the MPI version of the code is,

	/path/to/your/mpirun -np 100 program [OPTIONS] [FILE]  

Output
------
The program generates plain text files with names 'res_*.txt' to write 
simulation results. It generates one such file for each physical 
observables marked for measurement in the input. The full file 
name is derived from the observable name.

Adding a new model or lattice
-----------------------------

To add a new model or lattice, you will have to do some C++ coding.
Look at the files 'modellibrary.cc' and 'latticelibrary.cc'. 


Bug reports 
-----------
Send bug reports and suggestions to <amedhi@iisertvm.ac.in>.

Version history
---------------
cmc++ 1.0.0 - Sep 2016.

Copyright 
---------
Copyright (C) 2016-2016 by Amal Medhi, <amedhi@iisertvm.ac.in>.

All rights reserved.
