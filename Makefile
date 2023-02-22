#

include make.macro
EXEC= gebco_tool.exe gebco_xtrac.exe bathy_correction.exe

all: $(EXEC)


gebco_tool.exe : gebco_tool.f90
	$(F90)  gebco_tool.f90 -o gebco_tool.exe $(FFLAGS) 

gebco_xtrac.exe : gebco_xtrac.f90
	$(F90)  gebco_xtrac.f90 -o gebco_xtrac.exe $(FFLAGS) 

bathy_correction.exe : bathy_correction.f90
	$(F90)  bathy_correction.f90 -o bathy_correction.exe $(FFLAGS) 

clean:
	rm -f $(EXEC) *~


