
# For each program there must be a $(program_files) variable defined.
# The variable $(program_libs) should contain any other linking needs.
# The files variable defines the code needed to build
# the program.
Q=@
builddir=../build
srcdir=../srcs

# TODO Move these to includes.
GSLINCLUDE=-I/cluster/gsl/1.5/include
GSLLIBPATH=-L/cluster/gsl/1.5/lib
GSLLIBS=-lgsl -lgslcblas
GSLLINK=$(GSLLIBPATH) $(GSLLIBS)
MKLINCLUDE=-I/cluster/intel/mkl/10.0.1.014/include
MKLLIBPATH=-L/cluster/intel/mkl/10.0.1.014/lib/em64t
MKLLIBS=-lmkl_em64t -liomp5 -lpthread
MKLLINK=$(MKLLIBPATH) $(MKLLIBS)
LeapingMonkeyLink=-L../lib -lLeapingMonkey

COMMONLIBS=$(GSLLINK)
INCLUDES=$(GSLINCLUDE) $(MKLINCLUDE)

CC=icpc

PROGRAMS := TestSpectralHorizon TestLegendre \
	 TestConformal TestChebyshevExtrema TestChebyshevRoots\
	 TestBasis2D TestSphericalBasis TestExtrinsicCurvature\
	 TestHamiltonian TestTrumpet TestSingularPart TestFields\
	 TestNewtonRaphson1D TestSphericalHorizon

TestSphericalHorizon_files := TestSphericalHorizon.C SphericalHorizon.C\
	Fields.C Hamiltonian.C ExtrinsicCurvature.C SingularPart.C
TestSphericalHorizon_libs := $(LeapingMonkeyLink) $(MKLLINK)

TestNewtonRaphson1D_files := TestNewtonRaphson1D.C
TestNewtonRaphson1D_libs := $(LeapingMonkeyLink) $(MKLLINK)

TestSingularPart_files := TestSingularPart.C SingularPart.C ChebyshevRoots.C\
		Chebyshev.C Basis.C ExtrinsicCurvature.C
TestSingularPart_libs := $(MKLLINK)

TestTrumpet_files := TestTrumpet.C TrumpetConformal.C

TestSpectralHorizon_files := TestSpectralHorizon.C SpectralHorizon.C ConformalFactor.C
TestSpectralHorizon_libs := $(LeapingMonkeyLink) $(MKLLINK)

TestLegendre_files := TestLegendreMain.C TestLegendre.C Legendre.C Basis.C
TestLegendre_libs := $(MKLLINK)

TestConformal_files := TestConformal.C ConformalFactor.C

TestChebyshevExtrema_files := TestChebyshevExtrema.C ChebyshevExtrema.C\
			      Basis.C Chebyshev.C
TestChebyshevExtrema_libs := $(MKLLINK)

TestChebyshevRoots_files := TestChebyshevRoots.C ChebyshevRoots.C\
			      Basis.C Chebyshev.C
TestChebyshevRoots_libs := $(MKLLINK)

TestBasis2D_files := TestBasis2D.C Basis2D.C Basis.C ChebyshevRoots.C\
		     Chebyshev.C Legendre.C
TestBasis2D_libs := $(MKLLINK)

TestSphericalBasis_files := SphericalBasis.C TestSphericalBasis.C Basis.C\
			Chebyshev.C ChebyshevRoots.C Legendre.C Basis2D.C
TestSphericalBasis_libs := $(MKLLINK)

TestExtrinsicCurvature_files := TestExtrinsicCurvature.C ExtrinsicCurvature.C

TestHamiltonian_files := SphericalBasis.C Basis.C Chebyshev.C ChebyshevRoots.C\
			Basis2D.C Hamiltonian.C TestHamiltonian.C\
		       	ExtrinsicCurvature.C SingularPart.C
TestHamiltonian_libs := $(MKLLINK)

TestFields_files := TestFields.C Fields.C Hamiltonian.C ExtrinsicCurvature.C\
		SingularPart.C
TestFields_libs := $(LeapingMonkeyLink) $(MKLLINK)

TrumpetSolver_files := TrumpetSolver.C TrumpetConformal.C
TrumpetSolver_libs := $(LeapingMonkeyLink) $(MKLLINK)

MakeHorizons_files :=  MakeHorizons.C SpectralHorizon.C ConformalFactor.C
MakeHorizons_libs := $(LeapingMonkeyLink) $(MKLLINK)

Nonlinear_files := Nonlinear.C ExtrinsicCurvature.C Hamiltonian.C\
		SingularPart.C
Nonlinear_libs := $(LeapingMonkeyLink) $(MKLLINK)

#######################################################################
#
# Begin generic part. 
#
# Do not modify below here.
#
#######################################################################
.PHONY: clean all docs
all: $(LIBRARIES) $(PROGRAMS)

clean:
	@echo "----> clean"
	$(Q)rm -rf $(builddir)
	$(Q)rm -f $(PROGRAMS)

docs:
	$(Q)doxygen

# Functions to get the sources and objects from the individual
# filenames.
sources=$(patsubst %.C, $(srcdir)/%.C, $(1))
objs=$(patsubst %.C, $(builddir)/%.o, $(1))
link=$(COMMONLIBS) $$($(1)_libs)
prefix=.

libdir=$(prefix)/lib
$(libdir) :
	$(Q)mkdir $(libdir)

$(builddir):
	$(Q)mkdir $(builddir)

# This function takes the name of the program, and builds a Makefile
# rule to build the program. The objects are determined from the name of
# the program from the $(Program_files) variable. I more or less
# lifted this from the GNU make manual.
define PROGRAM_template
$(1): $$(call objs, $$($(1)_files))
	$(Q)$$(CC) $$^ $(call link,$(1)) -o $$@
	@echo "----> Built $(1)."
ALL_OBJS += $$(call objs, $$($(1)_files))
endef
define LIBRARY_template
$(1): $(libdir)/lib$(1).a

$(libdir)/lib$(1).a: $$(call objs, $$($(1)_files)) | $(libdir)
	$(Q)ar rcs $$@ $$^
	@echo "----> Built library: $(1)."
ALL_OBJS += $$(call objs, $$($(1)_files))
endef

# For each program, read in the Makefile rule to build it from the
# template.
$(foreach prog, $(PROGRAMS), $(eval $(call PROGRAM_template,$(prog))))
$(foreach lib, $(LIBRARIES), $(eval $(call LIBRARY_template,$(lib))))

# A single generic rule. Makes the .o and .d file from a .C file.
$(builddir)/%.o: $(srcdir)/%.C | $(builddir)
	$(Q)$(CC) $(INCLUDES) $(CPPFLAGS) -c $< -o $@
	$(Q)$(CC) $(INCLUDES) $(CPPFLAGS) -MT $@ -MM $< > $(builddir)/$*.d

# Include the dependency files, if they exist.
-include $(ALL_OBJS:.o=.d)

