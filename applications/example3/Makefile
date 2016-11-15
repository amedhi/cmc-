# Makefile options
include ../../root_dir.mk # must be included first
include ../../options.mk

CXX = $(CMC_CXX)
CXXFLAGS = $(CMC_CXXFLAGS)
LDFLAGS = $(CMC_LDFLAGS)
LIBS = $(CMC_LIBS)

#-------------------------------------------------------------
# Source files
SRCS = mc_ni2mnx.cc main.cc
HDRS = mc_ni2mnx.h
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
mc_ni2mnx.o : mc_ni2mnx.cc  mc_ni2mnx.h 
main.o : main.cc  mc_ni2mnx.h 

#-------------------------------------------------------
.PHONY: clean
clean:
	rm -f $(TAGT) *.o 
	
.PHONY: xclean
xclean:
	rm -f jobinfo* res_*.txt



