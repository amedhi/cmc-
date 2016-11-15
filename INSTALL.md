About
-------
The code provides the 'cmc++' library and sample classical Monte Carlo
applications developed using the library.  


Version
-------
cmc++ 1.0.0 - Oct 2016.

Build
------------
To build the library, first make a copy of the file 'make_options.mk' in the 
project root directory and name it as 'options.mk'. Edit this file to set your 
compilation options. Then do 'make' in the root directory. If successful, 
it builds the static library 'libcmc++.a' and installs it in the 'lib' directory
and the header files in the 'include' directory under the project root directory. 
There are two dependencies - Boost and Eigen C++. You will need latest C++ 
compilers as the code contains many C++11 features.  It it tested with 
latest g++ and clang++. 

To build the applications, go the the particular example directory under the
'applications' directory and do 'make' there.

Usage
-----
To run the application (say, ./a.out), do 

	./a.out [OPTIONS] [FILE]  

The program accepts a few optional arguments OPTIONS. Type 

	./a.out -help

to see the list of available options. The FILE argument, if present, 
is taken as the input filename, else a default filename 'input.parm' is 
looked for in the current directory.  
