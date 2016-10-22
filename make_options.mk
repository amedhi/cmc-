# User Configurable Options

# Building the cmc++ lib
#-------------------------------------------------------------
# 1. Set compiler option
#CMC_CXX=clang++ -std=c++11 # Clang compiler 
CMC_CXX=g++ -std=c++11 # GNU GCC compiler

ifndef CMC_CXX
$(error Makefile variable CMC_CXX not defined in options.mk; please define it.)
endif
#-------------------------------------------------------------
# 2. Compile flag options
# Flags to give the compiler for "release mode"
CMC_OPTFLAGS=-Wall -O3
# Flags to give the compiler for "debug mode"
CMC_DEBUGFLAGS=-DDEBUG_MODE -g -Wall -pedantic

#-------------------------------------------------------------
# 3. Boost and Eigen library
# Flags to give the compiler for "release mode"
BOOST_INCLUDE=-I/usr/local/include
EIGEN_INCLUDE=-I/usr/local/include
#BOOST_LIBS=-lboost_system -lboost_filesystem
#BOOST_LDFLAGS=-L/usr/local/lib

INCLUDE = $(BOOST_INCLUDE)
ifneq ($(BOOST_INCLUDE), $(EIGEN_INCLUDE))
INCLUDE += $(EIGEN_INCLUDE)
endif

CMC_CXXBFLAGS= $(CMC_OPTFLAGS) $(INCLUDE)
#-------------------------------------------------------------
# 4. Where to put the 'cmc' library & the includes
PREFIX=$(PROJECT_ROOT)
BUILD_DIR=$(PREFIX)/build
CMC_LIBDIR=$(PREFIX)/lib
CMC_INCLUDE=$(PREFIX)/include
CMC_CXXFLAGS= $(CMC_OPTFLAGS) $(INCLUDE) -I$(CMC_INCLUDE)
CMC_LDFLAGS=$(BOOST_LDFLAGS) -L$(CMC_LIBDIR)
CMC_LIBS=$(BOOST_LIBS) -lcmc++
#-------------------------------------------------------------
