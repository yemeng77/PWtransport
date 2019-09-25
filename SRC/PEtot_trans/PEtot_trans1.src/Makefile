all: PEtot_trans1

#include makefile.bassi
include makefile.intel

.SUFFIXES: .f .f90
.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(FC) $(F90FLAGS) -c $<

MODULES=load_data.f \
fft_data.f \
data.f

SRCL0=mainMV.f

SRCL1=moment.f

SRCL=Hpsi_comp.f\
Hpsi_comp_AllBandBP.f\
Proj1_comp.f\
Proj2_comp.f\
CG_comp.f\
DIIS_comp.f\
CG_AllBand.f\
ChebFD_BP.f\
lanczos_esti.f\
orth_comp.f\
orth_comp_DIIS.f\
getpot2L.f\
getpot3L.f\
getpot4L.f\
getpot4_force.f\
getpot5L.f\
getpot5_force.f\
GGAPBE.f\
d3fft_comp.f\
d3fft_comp_block.f\
cfft.f\
cfftd.f\
diag_comp.f\
djacobi.f\
convert_SLvr.f\
d3fft_real2.f\
d3fft_real2L.f\
d3fft_real2L2.f\
fftprep_comp.f\
fftprep_real2.f\
fftprep_real2L.f\
fftprep_real2L2.f\
fwdcpfft_comp.f\
fwdcpfft_comp_block.f\
fwdcpfft2.f\
fwdcpfft2L.f\
fwdcpfft2L2.f\
invcpfft_comp.f\
invcpfft_comp_block.f\
invcpfft2.f\
invcpfft2L.f\
invcpfft2L2.f\
gen_G_comp.f\
gen_G2_real.f\
gen_G2L_real.f\
gen_G2L2_real.f\
global_maxi.f\
global_sumr.f\
global_sumc.f\
fmin.f\
heapsort.f\
input.f\
init_ug_BP.f\
gaussj.f\
UxcCA.f\
UxcCA2.f\
w_line.f\
getwmask.f\
getwmaskX.f\
getwmask_q.f\
add_rho_beta.f\
getVrhoL.f\
mch_pulay.f\
mch_kerk.f\
Thomas3.f\
getNLsign.f\
getwq.f\
atomMV.f\
Etotcalc.f\
rhoIO.f\
ugIO.f\
ugIOBP.f\
wqIO.f\
beta_psiIO.f\
get_ALI.f\
gen_Gstar_ind.f\
symmop.f\
symmopf.f\
symmcheck.f\
forcLC.f\
forcNLq.f\
forcNLr.f\
ran1.f\
dens_out.f\
densWr_out.f\
readusp_head.f\
w_line_vwr.f\
w_line_usp.f\
LegendreSP.f\
clebsch_gordan.f\
get_Dij.f\
getwmask_dq.f\
get_VdqdR.f\
convert_2LtoL.f\
getV_Hartree.f\
getewald3D.f\
forcLC2.f\
getEextV.f\
getrho_only.f\
write_wg_BP.f\
occup.f


SRC_CG=orthogonal_symm.f90\
orthogonal_choleskyBP.f90\
orthogonal_choleskyBP2.f90\
orthogonal_projectionBP.f90\
orthogonal_projectionBP2.f90\
rotate_wfBP.f90\
rotate_wfBP_test.f90\
dot_product_BP.f90

OBJL=$(SRCL:.f=.o) $(SRC_obj) $(SRC_CG:.f90=.o)

OBJL0=$(SRCL0:.f=.o)

OBJL1=$(SRCL1:.f=.o)

OBJMOD=$(MODULES:.f=.o)

$(OBJL):param.escan_real $(MODULES)

$(OBJL0):param.escan_real $(MODULES)

$(OBJL1):param.escan_real $(MODULES)


PEtot_trans1: $(OBJMOD) $(OBJL0) $(OBJL) $(MODULES)
	$(LDR) -o $(EXEC_PETOT) $(OBJL0) $(OBJMOD) $(OBJL) $(LIB)

moment: $(OBJMOD) $(OBJL1) $(OBJL) $(MODULES) moment.f
	$(LDR) -o $(EXEC_PMOMENT) $(OBJL1) $(OBJMOD) $(OBJL) $(LIB)


clean:
	rm -f *.o *.mod *~
