os = $(shell uname -s)

INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(STARPICOPATH)
INCFLAGS      += -I./src


ifeq ($(os),Linux)
CXXFLAGS      = -O2 -fPIC -pipe -Wall -std=c++1z
CXXFLAGS     += -Wno-unused-variable
CXXFLAGS     += -Wno-unused-but-set-variable
CXXFLAGS     += -Wno-sign-compare
# # for gprof -- cannot combine with -g!
# CXXFLAGS     += -pg
# # for valgrind, gdb
# CXXFLAGS     += -g
else
CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
CXXFLAGS     += -Wno-return-type-c-linkage
## for debugging:
# CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       =
# # for gprof -- cannot combine with -g!
# LDFLAGS      += -pg
# # for valgrind, gdb
# LDFLAGS      += -g

LDFLAGSS      = --shared 
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = clang
endif

LDFLAGS	     += -lEG


# # uncomment for debug info in the library
# CXXFLAGS     += -g


ROOTLIBS      = $(shell root-config --libs)

LIBPATH       = $(ROOTLIBS) -L$(FASTJETDIR)/lib -L$(STARPICOPATH)
LIBS         += -lfastjet -lfastjettools -lTStarJetPico -lRecursiveTools

# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################
INCS = $(SDIR)/JetAnalyzer.hh
INCS = $(SDIR)/ppParameters.hh $(SDIR)/ppAnalysis.hh

###############################################################################
# standard rules
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $^ $(LIBS) -o $@

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all    : $(BDIR)/RunppAna
#	 doxy

$(SDIR)/dict.cxx 		: $(SDIR)/ktTrackEff.hh
	cd $(SDIR); rootcint -f dict.cxx -c -I. ./ktTrackEff.hh

$(ODIR)/dict.o 		: $(SDIR)/dict.cxx
$(ODIR)/ktTrackEff.o 	: $(SDIR)/ktTrackEff.cxx $(SDIR)/ktTrackEff.hh
$(ODIR)/JetAnalyzer.o   : ${SDIR}/JetAnalyzer.cxx ${INCS} ${SDIR}/JetAnalyzer.hh
$(ODIR)/ppAnalysis.o : $(SDIR)/ppAnalysis.cxx $(INCS) $(SDIR)/ppAnalysis.hh

# bin
$(BDIR)/RunppAna	:		$(ODIR)/RunppAna.o	$(ODIR)/JetAnalyzer.o	$(ODIR)/ppAnalysis.o  
###############################################################################
##################################### MISC ####################################
###############################################################################


doxy: html/index.html

html/index.html : $(INCS) src/* Doxyfile
#	doxygen
	@echo 
	@echo Updating documentation
	( cat Doxyfile ; echo "QUIET=YES" ) | doxygen -

clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -rvf $(BDIR)/*dSYM
	rm -rvf lib/*dSYM	
	rm -vf $(BDIR)/*
	rm -vf lib/*
	rm -vf $(SDIR)/dict.cxx $(SDIR)/dict.h

.PHONY : clean doxy
