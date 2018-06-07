
FC=ftn
LDR=ftn
FFLAGS =  -O3 -qfixed -qarch=pwr3 -qstrict -qspillsize=32648
LIB = $$SCALAPACK 

SRCL0=analyWave.f wave_decomp.f wave_system.f wave_electrode_3kpts.f wave_electrode.f 

OBJL0=$(SRCL0:.f=.o)

.f.o:
	$(FC)  -c $<

analyWave.x: $(OBJL0)  
	$(LDR)  -o analyWave.x $(OBJL0) $(LIB) 

clean:
	rm -f *.o
