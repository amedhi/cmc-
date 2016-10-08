#-------------------------------------------------------------
# Makefile for cmc++ library
#-------------------------------------------------------------
.PHONY: all
all: root_dir
	@cd src && $(MAKE)

.PHONY: root_dir
root_dir:
	@echo PROJECT_ROOT=`pwd` > root_dir.mk

.PHONY: install
install: all
	@cd src && $(MAKE) install

.PHONY: clean
clean:
	@cd src && $(MAKE) clean

.PHONY: bclean
bclean:
	@cd src && $(MAKE) bclean

.PHONY: distclean
distclean: clean
	@rm -f root_dir.mk 


