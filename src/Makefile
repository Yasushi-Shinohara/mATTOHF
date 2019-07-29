include ./make.inc

PROG = 1D-TDHF
OBJS = Constants.o Global_variables.o Functions.o \
	Initialization.o Prep_periodic_ion_potential.o \
	GS_calculation.o RT_calculation.o 1D-TDHF.o

$(PROG): $(OBJS)
	$(FC) -o $(PROG) $(OBJS) $(LN)

Constants.o:Constants.f90
	$(FC) -c Constants.f90 $(LN)
Global_variables.o:Global_variables.f90
	$(FC) -c Global_variables.f90 $(LN)
Functions.o:Functions.f90
	$(FC) -c Functions.f90 $(LN)
Initialization.o:Initialization.f90
	$(FC) -c Initialization.f90 $(LN)
Prep_periodic_ion_potential.o:Prep_periodic_ion_potential.f90
	$(FC) -c Prep_periodic_ion_potential.f90 $(LN)
GS_calculation.o:GS_calculation.f90
	$(FC) -c GS_calculation.f90 $(LN)
RT_calculation.o:RT_calculation.f90
	$(FC) -c RT_calculation.f90 $(LN)
1D-TDHF.o:1D-TDHF.f90
	$(FC) -c 1D-TDHF.f90 $(LN)
#
#Targets
#-----------------------------------------------------
clean:
	rm  1D-TDHF -f *~  *.o  *.mod
