# Makefile for the DMRG project
#-------------------------------------------------------------
CC = clang++
#CC = /opt/local/bin/g++
#CC = /opt/local/bin/g++-mp-5

#CPPFLAGS =-DWITH_LAPACK
OPTFLAGS = -O3 -Wall
#OPTFLAGS  = -g -O0
DBGFLAGS = -g -O0 -DDEBUG_MODE
MKLPATH = /opt/intel/mkl/lib
MKLINCLUDE = /opt/intel/mkl/include/intel64/lp64
#CCINCLUDE = $(PREFIX)/include
INCLUDE  = -I/usr/local/include -I $(MKLINCLUDE) 
CCFLAGS = -std=c++11 $(OPTFLAGS) $(INCLUDE) $(CPPFLAGS) 
CCGFLAGS = $(DBGFLAGS) $(INCLUDE) $(CPPFLAGS)  
LDFLAGS = $(OPTFLAGS) -L/usr/local/lib  -L$(MKLPATH) 
LDGFLAGS = $(DBGFLAGS) -L/usr/local/lib -L$(MKLPATH) #-L$(LIBGALAHAD)   
LIBS = -lboost_system -lboost_filesystem -lmkl_intel_lp64 -lmkl_sequential \
	-lmkl_core -lmkl_blas95_lp64 -lmkl_lapack95_lp64 

#-------------------------------------------------------------
TAGT = a.out
GTAGT = dbg.out
SRCS = cmdargs.cc 
SRCS+= inputparams.cc 
SRCS+= taskparams.cc 
SRCS+= worker.cc 
SRCS+= master_scheduler.cc
SRCS+= scheduler.cc
#-----------------------
SRCS+= expression.cc tokens.cc functions.cc objects.cc
#-----------------------
SRCS+= lattice.cc
SRCS+= latticelibrary.cc
SRCS+= graph.cc
SRCS+= qn.cc
SRCS+= quantum_operator.cc
SRCS+= sitebasis.cc
SRCS+= hamiltonian_term.cc
SRCS+= model.cc
SRCS+= modellibrary.cc
SRCS+= random.cc
SRCS+= observable_operator.cc
SRCS+= observables.cc
SRCS+= simulator.cc
SRCS+= measurement.cc
SRCS+= main.cc

HDRS = optionparser.h cmdargs.h inputparams.h worker.h task.h scheduler.h \
       expression.h shunting_yard.h tokens.h functions.h objects.h pack.h \
       constants.h lattice.h graph.h \
       qn.h sitebasis.h modelparams.h hamiltonian_term.h model.h \
       mcdata.h \
       random.h sitebasisstate.h observable_operator.h observables.h \
       simulator.h 

AUX = Makefile input.parm changelog
#-------------------------------------------------------------
# compilation and linking
VPATH = ./:./scheduler:./expression:./lattice:./model:./mcdata:./montecarlo

OBJS = $(patsubst %.cc,%.o, $(SRCS))
GOBJS= $(patsubst %.cc,debug_objs/%.o, $(SRCS))

all: $(TAGT)
$(TAGT) : $(OBJS)
	$(CC) -o $(TAGT) $(OBJS) $(LDFLAGS) $(LIBS)

$(GTAGT) : $(GOBJS)
	$(CC) -o $(GTAGT) $(GOBJS) $(LDFLAGS) $(LIBS)

%.o: %.cc
	$(CC) -c $(CCFLAGS) -o $@ $<

debug_objs/%.o: %.cc mkdebugdir
	$(CC) -c $(CCGFLAGS) -o $@ $<

mkdebugdir:
	mkdir -p debug_objs
	
#-------------------------------------------------------
# Detailed dependencies
DEPHDRS = optionparser.h cmdargs.h
cmdargs.o: $(DEPHDRS)
debug_objs/cmdargs.o: $(DEPHDRS)
DEPHDRS += inputparams.h
inputparams.o: inputparams.h
debug_objs/inputparams.o: $(DEPHDRS)
taskparams.o: inputparams.h
debug_objs/taskparams.o: $(DEPHDRS)
DEPHDRS += worker.h
DEPHDRS += task.h
worker.o: $(DEPHDRS)
debug_objs/task.o: $(DEPHDRS)
DEPHDRS += scheduler.h 
master_scheduler.o: $(DEPHDRS) 
scheduler.o: $(DEPHDRS) 
debug_objs/scheduler.o: $(DEPHDRS)
DEPHDRS += scheduler.h 
DEPHDRS += constants.h 
DEPHDRS += expression.h shunting_yard.h tokens.h functions.h objects.h 
expression.o : expression.h tokens.h functions.h objects.h 
DEPHDRS += lattice.h 
DEPHDRS += graph.h 
lattice.o: lattice.h 
latticelibrary.o: lattice.h 
graph.o: lattice.h graph.h
#DEPHDRS += blochbasis.h 
#blochbasis.o: blochbasis.h 
#model.o: $(DEPHDRS)
#modellibrary.o: $(DEPHDRS)
#DEPHDRS += hamiltonian.h 
#hamiltonian.o: $(DEPHDRS)
DEPHDRS += qn.h 
DEPHDRS += quantum_operator.h 
quantum_operator.o : quantum_operator.h 
qn.o: qn.h
DEPHDRS += sitebasis.h 
sitebasis.o: $(DEPHDRS)
DEPHDRS += modelparams.h hamiltonian_term.h 
hamiltonian_term.o: $(DEPHDRS)
DEPHDRS += model.h 
model.o: $(DEPHDRS)
modellibrary.o: $(DEPHDRS)
DEPHDRS += random.h 
random.o: random.h
DEPHDRS += sitebasisstate.h 
DEPHDRS += observable_operator.h 
observable_operator.o: observable_operator.h 
DEPHDRS += mcdata.h 
DEPHDRS += observables.h 
observables.o: $(DEPHDRS)
DEPHDRS += simulator.h 
simulator.o: $(DEPHDRS)
measurement.o: $(DEPHDRS)
main.o: $(DEPHDRS)

#-------------------------------------------------------
.PHONY: clean
clean:
	rm -f *.o debug_objs/*.o 
	
.PHONY: xclean
xclean:
	rm -f jobinfo* *.data*

