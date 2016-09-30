Version
-------
cmc++ 1.0.0 - Sep 2016.

Build
------------

There are two dependencies - Boost and Eigen C++. You will also need 
latest C++ compilers as the code contains many C++11 features. It works with
latest g++ and clang++. It should work with others also. See at the Makefile.

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
