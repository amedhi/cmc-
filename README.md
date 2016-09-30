Documentation                                               {#mainpage}
=============

Version
-------
cmc++ 1.0.0 - Sep 2016.

What is new
-----------
Start version ready.

Introduction
------------

cmc++ is a program written in C++ for doing Monte Carlo simulation for 
'classical' model Hamiltonians like the Ising model. 
The core codes are model and lattice geometry independent. This means 
it should work for all 'classical' models and all lattice geometries once 
these definitions are incorporated into the project in an appropriate manner. 
#The program can be build in serial mode or in MPI mode by simply a choice of
#compiler directive. In parallel mode, the program creates multiple 
#independent Markov chains in each processor and the final results 
#averaged over the results from these independent chains.  
The code is designed 
to make it easy to add new features - such as, adding a new model or lattice, 
adding a new physical observable to calculate or a new input parameter.

Compilation
------------

There are two dependencies - Boost and Eigen C++. You will also need 
latest C++ compilers as the code contains many C++11 features. It works with
latest g++ and clang++. It should work with others also.

Usage
-----
Once the executable is built and installed, run it as 

	program [OPTIONS] [FILE]  

Here 'program' is the name of the executable.  The program accepts a few 
optional arguments OPTIONS. Type 

	program -help

to see the list of available options. The FILE argument, if present, 
is taken as the input filename, else a default filename 'input.parm' is 
looked for in the current directory.  

Version history
---------------
cmc++ 1.0.0 - Sep 2016.

Copyright 
---------
Copyright (C) 2016-2016 by Amal Medhi, <amedhi@iisertvm.ac.in>.

All rights reserved.
