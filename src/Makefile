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
#If the files are not comiled using the genieric rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------

SRC := \
				main.f95 	      	\
				globals.f95	\
				string_module.f95	\
				errors_module.f95	\
				precision_module.f95	\
				input_module.f95	\
				output_module.f95

OBJ := $(SRC:.f95=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f95=.o)
MOD_FILTER := $(FILTER:.f95=.mod)

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
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f95
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f95
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f95] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
errors_module.o errors_module.mod: errors_module.f95 globals.mod
globals.o globals.mod: globals.f95
input_module.o input_module.mod: input_module.f95 globals.mod \
 errors_module.mod string_module.mod
main.o: main.f95 globals.mod input_module.mod output_module.mod
output_module.o output_module.mod: output_module.f95 globals.mod
precision_module.o precision_module.mod: precision_module.f95
string_module.o string_module.mod: string_module.f95 precision_module.mod