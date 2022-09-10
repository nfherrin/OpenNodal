# Defining variables
FC:=gfortran
EXE:=OpenNodal
EXT:=.exe

all :     TYPE :=
debug:    TYPE := _Debug
catch:    TYPE := _Catch
coverage: TYPE := _Coverage

all:      OPT := -O3 -Wall -Werror
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow -fcheck=all
catch:    OPT := -O0 -g -Wall -W -Wsurprising -Wextra #-Werror
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not compiled using the generic rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------

SRC := main.f90 \
	   edits_module.f90 \
	   globals.f90 \
	   string_module.f90 \
	   errors_module.f90 \
	   precisions.f90 \
	   xs_types.f90 \
	   input_module.f90 \
	   output_module.f90 \
	   solvers_module.f90 \
	   linalg_module.f90 \
	   analytic_module.f90

OBJ := $(SRC:.f90=.o)

#lapack libraries
BLIBS = -llapack -lblas

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all debug catch coverage clean reset

EXEC_NAME=$(EXE)
EXEC_FILE=../$(EXEC_NAME)

all: $(EXEC_FILE)

debug: $(EXEC_FILE)

catch: $(EXEC_FILE)

coverage: $(EXEC_FILE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	$(RM) $(OBJ) *.mod

#Intended to reset directory to fresh state with no exe or artifacts
reset : clean
	$(RM) ../$(EXE)*

#Generics for source files
%.o : %.f90
	$(FC) $(OPT) -c $<

#Executable rule
$(EXEC_FILE): $(OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $^ $(BLIBS)

#-----------------------------------------------------------------------------------
#Dependency List
#-----------------------------------------------------------------------------------
analytic_module.o : errors_module.o globals.o precisions.o string_module.o xs_types.o
edits_module.o : globals.o string_module.o xs_types.o
errors_module.o : globals.o
globals.o : precisions.o
input_module.o : edits_module.o errors_module.o globals.o string_module.o
linalg_module.o : precisions.o
main.o : analytic_module.o errors_module.o globals.o input_module.o output_module.o solvers_module.o string_module.o xs_types.o
output_module.o : errors_module.o globals.o string_module.o xs_types.o
solvers_module.o : errors_module.o globals.o string_module.o xs_types.o
string_module.o : precisions.o
xs_types.o : precisions.o
