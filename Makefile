# TODO Move these to includes.
GSLINCLUDE=-I/cluster/gsl/1.5/include
GSLLIBPATH=-L/cluster/gsl/1.5/lib
GSLLIBS=-lgsl -lgslcblas
GSLLINK=$(GSLLIBPATH) $(GSLLIBS)
MKLINCLUDE=-I/cluster/intel/mkl/10.0.1.014/include
MKLLIBPATH=-L/cluster/intel/mkl/10.0.1.014/lib/em64t
MKLLIBS=-lmkl_em64t -liomp5 -lpthread
MKLLINK=$(MKLLIBPATH) $(MKLLIBS)
LeapingMonkeyLink=-Llib -lLeapingMonkey

CC=icpc

# For each program there must be a $(program_files) variable defined.
# The variable $(program_libs) should contain any other linking needs.
# The files variable defines the code needed to build
# the program.
Q=@
builddir=build
srcdir=srcs
COMMONLIBS=$(GSLLINK)
INCLUDES=$(GSLINCLUDE) $(MKLINCLUDE)
PROGRAMS := MakeHorizons Nonlinear MakeFields
LIBRARIES := LeapingMonkey

MakeHorizons_files :=  MakeHorizons.C SpectralHorizon.C ConformalFactor.C
MakeHorizons_libs := $(LeapingMonkeyLink) $(MKLLINK)

Nonlinear_files := Nonlinear.C ExtrinsicCurvature.C Hamiltonian.C\
		SingularPart.C
Nonlinear_libs := $(LeapingMonkeyLink) $(MKLLINK)

MakeFields_files := MakeFields.C FieldMaker.C Fields.C Hamiltonian.C\
		ExtrinsicCurvature.C SingularPart.C SphericalHorizon.C\
		TrumpetConformal.C
MakeFields_libs := $(LeapingMonkeyLink) $(MKLLINK)

LeapingMonkey_files := Basis.C ChebyshevRoots.C SphericalBasis.C \
                Basis2D.C Chebyshev.C Legendre.C NewtonRaphson1D.C \
		ChebyshevExtrema.C

#######################################################################
#
# Begin generic part. 
#
# Do not modify below here.
#
#######################################################################
.PHONY: clean all docs install
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
depfiles=$(patsubst %.C, $(builddir)/%.d, $(1))
link=$(COMMONLIBS) $$($(1)_libs)
prefix=.

libdir=$(prefix)/lib
incdir=$(prefix)/include
$(libdir) :
	$(Q)mkdir $(libdir)

$(incdir) :
	$(Q)mkdir $(incdir)

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

$(libdir)/lib$(1).a: $$(call objs, $$($(1)_files)) | $(libdir) $(incdir)
	$(Q)ar rcs $$@ $$^
	$(Q)cp `cat $$(call depfiles, $$($(1)_files)) | sed 's/ /\n/g' | grep $(srcdir) | grep '\.h' | sort | uniq | sed 's/\n/ /g'` $(incdir)
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

