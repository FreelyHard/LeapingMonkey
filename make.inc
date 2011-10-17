LAPACKINCLUDE=-Ie:\include
LAPACKLIBPATH=-Le:\lib
COMMONLIBS=
LAPACKLIBS=-llapacke -llapack_WIN7 -lcblas_WIN7 -lblas_WIN7 -lgfortran
LAPACKLINK=$(LAPACKLIBPATH) $(LAPACKLIBS)
LeapingMonkeyLink=-Le:\lib -lLeapingMonkey
INCLUDES=$(LAPACKINCLUDE)

# Operating system dependent stuff.
DIRSEP=\\
rm=del
cp=copy

CC=g++
LD=gfortran
