# Makefile options
include ../../root_dir.mk # must be included first
include ../../options.mk

CXX = $(CMC_CXX)
CXXFLAGS = $(CMC_CXXFLAGS)
LDFLAGS = $(CMC_LDFLAGS)
LIBS = $(CMC_LIBS)

#-------------------------------------------------------------
# Target
TAGT = nimnx
all: $(TAGT)
nimnx: mc_nimnx.o libcmc++.a
	$(CXX) -o nimnx mc_nimnx.o $(LDFLAGS) $(LIBS)

mc_nimnx.o: mc_nimnx.cc
	$(CXX) -c $(CXXFLAGS) -o mc_nimnx.o mc_nimnx.cc

libcmc++.a: 
	@cd ../../ && $(MAKE)

#-------------------------------------------------------
.PHONY: clean
clean:
	rm -f $(TAGT) *.o 
	
.PHONY: xclean
xclean:
	rm -f jobinfo* res_*.txt



