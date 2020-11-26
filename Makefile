SRC      = .

FF       = gfortran

#BCHECK   = -fbounds-check
OPTIM = -o2

FFLAGS   = $(OPTIM) $(BCHECK) 
LFLAGS   = $(OPTIM) 

OBJ      = diskevol.o diskevol_module.o dustmodule.o natconst_module.o nrecip_module.o 

all:	            diskevol

diskevol:           $(OBJ) Makefile
	            $(FF) $(LFLAGS) $(OBJ) $(LIBS) -o $@ 

diskevol.o:         diskevol.f90 diskevol_module.o dustmodule.o natconst_module.o Makefile
	            $(FF) -c $(FFLAGS) diskevol.f90 -o $@

diskevol_module.o:  $(SRC)/diskevol_module.f90 dustmodule.o nrecip_module.o natconst_module.o Makefile
	            $(FF) -c $(FFLAGS) $(SRC)/diskevol_module.f90 -o $@

dustmodule.o: dustmodule.f90 natconst_module.o Makefile
	            $(FF) -c $(FFLAGS) dustmodule.f90 -o $@

nrecip_module.o:    $(SRC)/nrecip_module.f90 Makefile
	            $(FF) -c $(FFLAGS) $(SRC)/nrecip_module.f90 -o $@

natconst_module.o:  $(SRC)/natconst_module.f90 Makefile
	            $(FF) -c $(FFLAGS) $(SRC)/natconst_module.f90 -o $@


clean:
	@rm -f	*.o *.mod 
	@echo OBJECT and MODULE files removed.

cleanall:
	@rm -f	*.o *.mod *~ *.dat *.info *.pyc fort.321 diskevol
	@echo Cleaned directory back to basic
