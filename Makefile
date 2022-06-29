##############################################################################
# D310 MAKEFILE
##############################################################################


OBJS=main.o\
input.o\
matrix3.o\
gfun.o\
cell.o\
cellFile.o\
cellPROFESS.o\
cellVASP.o\
cellQE.o\
cellXYZ.o\
cellLAMMPS.o\
cellABACUS.o\
cellRAW.o\
cellQE2.o\
cellPWmat.o\
atoms.o\
ext.o\
pdf.o\
pdf_added.o\
pdf5.o\
pdf2d.o\
ssf.o\
ssf_selected.o\
dsf.o\
vel.o\
velcor.o\
powers.o\
vacuum.o\
data3D.o\
math.o\
pseudo.o\
iprof.o\
isf.o\
isf2.o\
bdf.o\
bdf_rcut.o\
average.o\
insert.o\
ili.o\
ili_3D.o\
mdp.o\
mdp2.o\
mdp3.o\
HBs.o\
HB_stat.o\
mj.o \
tetra_order.o\
water.o\
movie.o\
wannier.o\
wannier1.o\
eig.o\
void.o\
hyper.o\
trajadj.o\
dist.o\
pre.o\
msd.o\
msd_multiple.o\
Honeycutt.o\
waterwire.o\
waterwire2.o\
ww_compress.o\
xy_profile.o\
tune_stru.o\
reorganize.o\
directional.o\
dielectric.o\
xsf.o\
density2D.o\
HBs_near_PT.o\
qe_input.o\
PT_snapshot.o\
PT_snapshot2.o\
fp_check.o\
band_gap.o\
incremental_pdf.o\
special_msd.o\
stress_average.o\
mass_center.o\
HB_angle.o\
OH_movie.o\
movie2.o\
HB_stat2.o\
Wan_centers_stat.o\
dist2.o\
incremental_pdf2.o\
incremental_pdf3.o\
HB_stat3.o\
HB_stat4.o\
AngularJump.o\
HBs_near_AngularJump.o\
HBs_near_DoubleDonor.o\
DoubleDonor.o\
HBs_near_TwistThird.o\
HB_correlation.o\
HB_correlation2.o\
nonHB_correlation.o\
nonHB_correlation2.o\
nonHB_correlation3.o\
angular_correlation.o\
HB_break.o\
oho_angle.o\
orientation_tcf.o\
bdf_rcut1.o\
first_shell_angle.o\
#ring_group.o\


#CC=icc -g
#CC=mpiicc 
#CC=mpicxx
#CC=CC
#HONG=-D__MPI
#CC=mpiicc
CC=g++
OPTION=-O3 

.cpp.o:
	${CC} ${OPTION} -c ${HONG} $< -o $@

candela :${OBJS}
	${CC} ${OPTION} -o candela.exe ${OBJS}

clean:
	rm -f *.o *.exe
