# Defining variables
FC:=gfortran
EXE:=OpenNodal
EXT:=.exe

all:			TYPE :=
debug:    TYPE := _Debug
catch:    TYPE := _Catch
coverage: TYPE := _Coverage

all:      OPT := -O3
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow
catch:    OPT := -O0 -g -Wall -W -Wsurprising -Werror
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not compiled using the generic rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------

SRC := \
				main.f90 	      	\
				edits_module.f90 \
				globals.f90	\
				string_module.f90	\
				errors_module.f90	\
				precisions.f90	\
				xs_types.f90	\
				input_module.f90	\
				output_module.f90	\
				solvers_module.f90

OBJ := $(SRC:.f90=.o)
MOD := $(OBJ:.o=.mod) # TODO this is wrong, not all f90 have a mod (e.g., main)

FILTER := \

OBJ_FILTER := $(FILTER:.f90=.o)
MOD_FILTER := $(FILTER:.f90=.mod)

#lapack libraries
BLIBS = -llapack -lblas

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all debug catch coverage clean reset

print-%  : ; @echo $* = $($*)

all: $(EXE)

debug: $(EXE)

catch: $(EXE)

coverage: $(EXE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	rm -f $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	rm -f $(MOD) $(CONT_MOD) $(LIB_MOD)
	rm -f $(COMP_DEP)

#Intended to reset directory to fresh state with no exe or artifacts
reset: clean
	rm -f *.gcno *.gcda *.o *.mod
	rm -f -r $(EXE)*.dSYM
	rm -f ../$(EXE)*

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f90
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f90
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) $(CONT_OBJ) $(LIB_OBJ) $(BLIBS)
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f90] repeatedly until clean compile to update rules below
# TODO makefile is busted.
# Can't use "mod" in dependency list. Timestamp doesn't update. Need to use "o"
#-----------------------------------------------------------------------------------
edits_module.o edits_module.mod: edits_module.f90 string_module.mod \
 xs_types.mod precisions.mod globals.mod errors_module.mod
errors_module.o errors_module.mod: errors_module.f90 globals.mod
globals.o globals.mod: globals.f90 precisions.mod xs_types.mod
input_module.o input_module.mod: input_module.f90 globals.mod \
 errors_module.mod string_module.mod edits_module.mod
main.o: main.f90 errors_module.mod globals.mod input_module.mod \
 output_module.mod solvers_module.mod
output_module.o output_module.mod: output_module.f90 globals.mod
precisions.o precisions.mod: precisions.f90
solvers_module.o solvers_module.mod: solvers_module.f90 globals.mod \
 errors_module.mod precisions.mod globals.mod precisions.mod \
 string_module.mod
string_module.o string_module.mod: string_module.f90 precisions.mod
xs_types.o xs_types.mod: xs_types.f90 precisions.mod