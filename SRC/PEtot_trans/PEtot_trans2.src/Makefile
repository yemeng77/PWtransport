all: PEtot_trans2

.SUFFIXES : .o .f90 .f
#include makefile.franklin
include makefile.intel

MODULES=load_data.f \
fft_data.f \
data.f 

SRCL0=mainMV.f

SRCL=Hpsi_comp.f\
CG_comp.f\
CG_new.f\
getpot2.f\
getpot3.f\
getpot4.f\
getpot4_force.f\
getpot5.f\
getpot5_force.f\
GGAPBE.f\
d3fft_comp.f\
cfft.f\
cfftd.f\
diag_comp.f\
djacobi.f\
d3fft_real2.f\
fftprep_comp.f\
fftprep_real2.f\
fwdcpfft_comp.f\
fwdcpfft2.f\
invcpfft_comp.f\
invcpfft2.f\
gen_G_comp.f\
gen_G2_real.f\
global_maxi.f\
global_sumr.f\
global_sumc.f\
fmin.f\
heapsort.f\
input.f\
init_ug.f\
gaussj.f\
UxcCA.f\
UxcCA2.f\
w_line.f\
getewald.f\
getEextV.f\
add_slab_dipV.f\
getwmask.f\
getwmaskX.f\
getVrho.f\
getrho_only.f\
get_atomW.f\
mch_pulay.f\
mch_kerk.f\
Thomas3.f\
getNLsign.f\
getwq.f\
atomMV.f\
Etotcalc.f\
Etranscalc.f\
CG_linear.f\
find_wrl.f\
wave_electrode_3kpts_LW.f\
rhoIO.f\
rhoIO_comp.f\
get_ALI.f\
ugIO.f\
wqIO.f\
occup.f\
gen_Gstar_ind.f\
symmop.f\
symmopf.f\
symmcheck.f\
forcLC.f\
forcNLq.f\
forcNLr.f\
ran1.f\
dens_out.f\
dens_outWG.f\
dens_outWGsp.f\
system_orth_comp.f\
system_orth_comp_N.f\
system_czheev.f\
system_flush.f\
system_time.f\
lanczos_esti.f\
QMR_linear.f

MOD_MINRESQLP=zminresqlpModule.f90

SRC_MINRESQLP=MINRESQLP_linear.f

OBJL=$(SRCL:.f=.o) $(SRC_obj)

OBJL0=$(SRCL0:.f=.o)

OBJ_MINRESQLP=$(SRC_MINRESQLP:.f=.o)

OBJMOD=$(MODULES:.f=.o)

OBJMOD_MINRESQLP=$(MOD_MINRESQLP:.f90=.o)

$(OBJL):param.escan_real $(OBJMOD)

$(OBJL0):param.escan_real $(OBJMOD)

$(OBJ_MINRESQLP):param.escan_real $(OBJMOD) $(OBJMOD_MINRESQLP)

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(FFLAGS) -c $<

PEtot_trans2: $(OBJL0) $(OBJMOD) $(OBJL) $(OBJ_MINRESQLP) $(OBJMOD_MINRESQLP)
	$(LDR)  -o PEtot_trans2 $^ $(LIB)

clean:
	rm -f *.o *.mod
