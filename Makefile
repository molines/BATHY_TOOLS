#

include make.macro
EXEC= gebco_tool.exe gebco_xtrac.exe

all: $(EXEC)


gebco_tool.exe : gebco_tool.f90
	$(F90)  gebco_tool.f90 -o gebco_tool.exe $(FFLAGS) 

gebco_xtrac.exe : gebco_xtrac.f90
	$(F90)  gebco_xtrac.f90 -o gebco_xtrac.exe $(FFLAGS) 

clean:
	rm -f $(EXEC) *~


