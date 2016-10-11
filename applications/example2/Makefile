# Makefile options
include ../../root_dir.mk # must be included first
include ../../options.mk

CXX = $(CMC_CXX)
CXXFLAGS = $(CMC_CXXFLAGS)
LDFLAGS = $(CMC_LDFLAGS)
LIBS = $(CMC_LIBS)

#-------------------------------------------------------------
# Source files
SRCS = mcsimulation.cc main.cc
HDRS = mcsimulation.h
#-------------------------------------------------------------
# Target
TAGT=a.out

OBJS = $(patsubst %.cc,%.o, $(SRCS))

all: $(TAGT)
$(TAGT) : $(OBJS)
	$(CXX) -o $(TAGT) $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

# detail dependencies
mcsimulation.o : mcsimulation.cc  mcsimulation.h 
main.o : main.cc  mcsimulation.h 

#-------------------------------------------------------
.PHONY: clean
clean:
	rm -f *.o 
	
.PHONY: xclean
xclean:
	rm -f $(TAGT) jobinfo* res_*.txt



