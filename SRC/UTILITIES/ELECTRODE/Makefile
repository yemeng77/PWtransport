#Makefile for connect(_mpi).x

FC=ftn
LDR=$(FC)
FFLAGS = -g -traceback

std: connect.x
mpi: connect_mpi.x
all: std mpi
.PHONY : std mpi all

.f.o:
	$(FC) $(FFLAGS) -c $<

connect.x: find_evan.o connect.o
	$(LDR) -o $@ $^
	
connect_mpi.x: find_evan.o connect_mpi.o
	$(LDR) -o $@ $^

clean:
	rm -f *.o

clean-all: clean
	rm -f *.x