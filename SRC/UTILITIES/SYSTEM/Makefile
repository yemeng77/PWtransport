
FC=ftn -g -traceback
LDR=ftn
LIB =  

SRCL0=find_wrl.f wave_electrode_3kpts_LW.f wave_electrode_interp.f

OBJL0=$(SRCL0:.f=.o)

.f.o:
	$(FC) -c $<

find_wrl.x: $(OBJL0)  
	$(LDR) -o find_wrl.x $(OBJL0) $(LIB) 

clean:
	rm -f *.o
