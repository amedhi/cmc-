#-------------------------------------------------------------
# Makefile for cmc++ library
#-------------------------------------------------------------
include ../root_dir.mk # must be included first
include ../options.mk
#-------------------------------------------------------------
# Source files
CMC_SRCS = scheduler/cmdargs.cc 
CMC_SRCS+= scheduler/inputparams.cc 
CMC_SRCS+= scheduler/taskparams.cc 
CMC_SRCS+= scheduler/worker.cc 
CMC_SRCS+= scheduler/master_scheduler.cc
CMC_SRCS+= scheduler/scheduler.cc
CMC_SRCS+= expression/expression.cc 
CMC_SRCS+= expression/tokens.cc 
CMC_SRCS+= expression/functions.cc
CMC_SRCS+= expression/objects.cc
CMC_SRCS+= lattice/lattice.cc
CMC_SRCS+= lattice/latticelibrary.cc
CMC_SRCS+= lattice/graph.cc
CMC_SRCS+= model/qn.cc
CMC_SRCS+= model/quantum_operator.cc
CMC_SRCS+= model/sitebasis.cc
CMC_SRCS+= model/hamiltonian_term.cc
CMC_SRCS+= model/model.cc
CMC_SRCS+= model/modellibrary.cc
CMC_SRCS+= observable/observables.cc
CMC_SRCS+= montecarlo/random.cc
CMC_SRCS+= montecarlo/observable_operator.cc
#CMC_SRCS+= montecarlo/measurement.cc
CMC_SRCS+= montecarlo/simulator.cc
#-------------------------------------------------------------
# Headers
CMC_HDRS=scheduler/optionparser.h scheduler/cmdargs.h \
         scheduler/inputparams.h scheduler/worker.h scheduler/task.h \
         scheduler/scheduler.h \
         expression/expression.h expression/shunting_yard.h \
         expression/tokens.h expression/functions.h expression/objects.h \
         expression/pack.h \
         lattice/constants.h lattice/lattice.h lattice/graph.h \
         model/qn.h model/quantum_operator.h model/sitebasis.h model/modelparams.h \
         model/hamiltonian_term.h model/model.h \
         observable/mcdata.h observable/observables.h \
         montecarlo/random.h montecarlo/sitebasisstate.h \
	 montecarlo/observable_operator.h montecarlo/simulator.h \
	 simulation.h
#-------------------------------------------------------------
# Target
TGT=$(CMC_LIBDIR)/libcmc++.a

# Put all auto generated stuff to this build dir.
ifeq ($(BUILD_DIR), $(CURDIR))
$(error In-source build is not allowed, choose another build directory)
endif

# All .o files go to BULD_DIR
OBJS=$(patsubst %.cc,$(BUILD_DIR)/%.o,$(CMC_SRCS))
# Gcc/Clang will create these .d files containing dependencies.
DEPS=$(patsubst %.o,%.d,$(OBJS)) 
# compiler flags
CXXFLAGS= $(CMC_INCLUDEFLAGS) $(OPTIMIZATIONS) 

CMC_INCLDIR=$(CMC_INCLUDE)/cmc
INCL_HDRS = $(addprefix $(CMC_INCLDIR)/, $(CMC_HDRS))

.PHONY: all
all: $(TGT) $(INCL_HDRS)

$(TGT): $(OBJS)
	@echo "Building static library $(TGT)" 
	@mkdir -p $(@D)
	@echo "ar cr $(@F) $(^F)"
	@ar cr $@ $(OBJS)
	ranlib $@

# Include all .d files
-include $(DEPS)

$(BUILD_DIR)/%.o: %.cc
	@mkdir -p $(@D)
	@echo "$(CMC_CXX) -c $(CMC_CXXBFLAGS) -o $(@F) $(<F)"
	@$(CMC_CXX) -MMD -c $(CMC_CXXBFLAGS) -o $@ $<

$(CMC_INCLDIR)/%.h: %.h 
	@mkdir -p $(@D)
	@echo "Copying $< to 'include'" 
	@cp -f $< $@

# installation
#prefix = ../install#/usr/local
#libdir = $(prefix)/lib
#includedir = $(prefix)/include/cmc++

.PHONY: install
install:	
	@echo "Already installed in $(CMC_LIBDIR) and $(CMC_INCLDIR)" 

.PHONY: clean
clean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
	@echo "Removing $(TGT) and headers in $(CMC_INCLDIR)/"
	@rm -f $(TGT) $(INCL_HDRS) 

.PHONY: bclean
bclean:	
	@echo "Removing temporary files in the build directory"
	@rm -f $(OBJS) $(DEPS) 
