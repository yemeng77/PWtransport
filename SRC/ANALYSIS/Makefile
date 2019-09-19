
FC=ftn -g -traceback
LDR=ftn
LIB =  

SRCL0=analyWave.f wave_decomp.f wave_system.f wave_electrode_3kpts.f wave_electrode.f wave_electrode_interp.f 

OBJL0=$(SRCL0:.f=.o)

SRCL1=analyWave_mpi.f find_state.f wave_decomp_ym.f wave_system_ym.f wave_electrode_3kpts.f wave_electrode.f wave_electrode_interp.f

OBJL1=$(SRCL1:.f=.o)

.f.o:
	$(FC) -c $<

analyWave.x: $(OBJL0)  
	$(LDR) -o analyWave.x $(OBJL0) $(LIB)

analyWave_mpi.x: $(OBJL1)  
	$(LDR) -o analyWave_mpi.x $(OBJL1) $(LIB)

mpi: analyWave_mpi.x

clean:
	rm -f *.o
