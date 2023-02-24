#

include make.macro
EXEC= gebco_tool.exe gebco_xtrac.exe bedmachine_mask.exe bathy_correction.exe bedmachine_idraft.exe

all: $(EXEC)


gebco_tool.exe : gebco_tool.f90
	$(F90)  gebco_tool.f90 -o gebco_tool.exe $(FFLAGS) 

gebco_xtrac.exe : gebco_xtrac.f90
	$(F90)  gebco_xtrac.f90 -o gebco_xtrac.exe $(FFLAGS) 

bedmachine_mask.exe : bedmachine_mask.f90
	$(F90)  bedmachine_mask.f90 -o bedmachine_mask.exe $(FFLAGS) 

bedmachine_idraft.exe : bedmachine_idraft.f90
	$(F90)  bedmachine_idraft.f90 -o bedmachine_idraft.exe $(FFLAGS) 

bathy_correction.exe : bathy_correction.f90
	$(F90)  bathy_correction.f90 -o bathy_correction.exe $(FFLAGS) 

clean:
	rm -f $(EXEC) *~


