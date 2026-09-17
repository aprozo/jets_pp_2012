# Builds bin/RunppAna, the Stage-1 jet-finding binary.
# Needs ROOTSYS, FASTJETDIR and STARPICOPATH set, so run it inside
# star_star.simg (scripts/runimage.sh, or `run_production.sh build`).

os = $(shell uname -s)
$(shell mkdir -p src/obj bin)

INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(STARPICOPATH)
INCFLAGS      += -I./src

ifeq ($(os),Linux)
CXXFLAGS      = -O2 -fPIC -pipe -Wall -std=c++1z
CXXFLAGS     += -Wno-unused-variable
CXXFLAGS     += -Wno-unused-but-set-variable
CXXFLAGS     += -Wno-sign-compare
else
CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
CXXFLAGS     += -Wno-return-type-c-linkage
endif

ifeq ($(os),Linux)
LDFLAGS       =
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
endif

ifeq ($(os),Linux)
CXX          = g++
else
CXX          = clang
endif

LDFLAGS	     += -lEG

ROOTLIBS      = $(shell root-config --libs)

LIBPATH       = $(ROOTLIBS) -L$(FASTJETDIR)/lib -L$(STARPICOPATH)
LIBS         += -lfastjet -lfastjettools -lTStarJetPico -lRecursiveTools

SDIR          = src
ODIR          = src/obj
BDIR          = bin

# Touching any of these headers rebuilds every object.
INCS = $(SDIR)/JetAnalyzer.hh $(SDIR)/ppParameters.hh $(SDIR)/ppAnalysis.hh $(SDIR)/JetQAHistogramManager.hh

$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $^ $(LIBS) -o $@

all    : $(BDIR)/RunppAna

$(ODIR)/JetAnalyzer.o   : ${SDIR}/JetAnalyzer.cxx ${INCS} ${SDIR}/JetAnalyzer.hh
$(ODIR)/JetQAHistogramManager.o : $(SDIR)/JetQAHistogramManager.cxx $(INCS) $(SDIR)/JetQAHistogramManager.hh
$(ODIR)/ppAnalysis.o : $(SDIR)/ppAnalysis.cxx $(INCS) $(SDIR)/ppAnalysis.hh

$(BDIR)/RunppAna	:		$(ODIR)/RunppAna.o	$(ODIR)/JetAnalyzer.o	$(ODIR)/JetQAHistogramManager.o	$(ODIR)/ppAnalysis.o

clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -rvf $(BDIR)/*dSYM
	rm -rvf lib/*dSYM	
	rm -vf $(BDIR)/*
	rm -vf lib/*.so lib/*.o

.PHONY : clean
