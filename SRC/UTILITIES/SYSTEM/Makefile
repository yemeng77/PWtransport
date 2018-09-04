
FC=ftn -g -traceback
LDR=ftn
LIB =  

SRCL0=find_state_mpi.f wave_electrode_3kpts_LW.f wave_electrode_interp.f write_wl.f

OBJL0=$(SRCL0:.f=.o)

.f.o:
	$(FC) -c $<

find_state.x: $(OBJL0)  
	$(LDR) -o find_state.x $(OBJL0) $(LIB) 

clean:
	rm -f *.o
