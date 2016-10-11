# Makefile options
include ../../root_dir.mk # must be included first
include ../../options.mk

CXX = $(CMC_CXX)
CXXFLAGS = $(CMC_CXXFLAGS)
LDFLAGS = $(CMC_LDFLAGS)
LIBS = $(CMC_LIBS)

#-------------------------------------------------------------
# Source files
SRCS = mcsimulation.cc 
#-------------------------------------------------------------
# Target
TAGT=a.out

OBJS = $(patsubst %.cc,%.o, $(SRCS))

all: $(TAGT)
$(TAGT) : $(OBJS)
	$(CXX) -o $(TAGT) $(OBJS) $(LDFLAGS) $(LIBS)

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) -o $@ $<

#-------------------------------------------------------
.PHONY: clean
clean:
	rm -f *.o 
	
.PHONY: xclean
xclean:
	rm -f $(TAGT) jobinfo* res_*.txt



