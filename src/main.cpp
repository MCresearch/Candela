#include "gfun.h"    // global functions
#include "input.h"   // input prameters
#include "write.h"	//write geo input file, etc.
#include "ext.h"     // extend the original cell to supercell
#include "vacuum.h"  // add vacuum to geometry
#include "vel.h"     // plot the distribution of velocties 
#include "pdf.h"     // pair distribution function
#include "pdf2d.h"     // pair distribution function on 2D plane
#include "pdf5.h"     // pair distribution function for the first few neighbours.
#include "ssf.h"     // static structure factor
#include "ssf_selected.h" // static structure factor calculated for selected vectors 
#include "dsf.h"     // dynamics structure factor
#include "ele_conductivity.h"   //qianrui for electric conductivity
#include "elecond_contribute.h" //calculate contribution of each band to conductivities
#include "velcor.h"  // velocity autocorrelation function
#include "powers.h"  // power spectra (13-03-28)
#include "data3D.h"  // arrange 3D data, such as density
#include "pseudo.h"  // generate pseudopotentials
#include "iprof.h"   // generate ionic density prfile
#include "isf.h"     // intermediate scattering functions
#include "isf2.h"		//qianrui new way to calculate intermediate scattering functions
#include "bdf.h"     // bond angle distribution functions
#include "bdf_rcut.h"  // bond angle distribution functions
#include "average.h" // average the ions and vels
#include "insert.h"  // insert atoms into existing geometry
#include "ili.h" // compute the instataneous liquid interfae and its gradient
#include "ili_3D.h" // compute the 3D data of ili
#include "movie.h" // generate a movie
#include "mdp.h" // compute the mean density profile based on instataneous liquid interface
#include "mdp2.h" // compute the mean density profile2 based on mean liquid interface
#include "mdp3.h" // compute the mean density profile3 based on mean liquid interface
#include "HBs.h" // compute the hydrogen bonds
#include "mj.h" // compute the hydrogen bonds
#include "tetra_order.h" // compute the tetrahedral order
#include "wannier.h" // compute the wannier function centers
#include "wannier1.h" // analyze the Wannier functions for CH systems
#include "eig.h" // compute the distribution of eigenvalues
#include "void.h" // create voids
#include "hyper.h" // hyper coordination
#include "trajadj.h" // trajectory-adjacent atoms 2D plot
#include "dist.h" // compute 3D atom distribution of a specfic species such as hydroxide
#include "pre.h" // compute the presolvation structure
#include "msd.h" // compute mean square displacement
#include "msd_multiple.h" // compute multiple mean square displacement
#include "Honeycutt.h" // compute HoneycuttAnderson Short-Range Order
#include "waterwire.h" // compute the waterwire
#include "waterwire2.h" // compute the waterwire
#include "ww_compress.h" // compute the waterwire compression
#include "xy_profile.h" // compute the location of ion from xy plane.
#include "tune_stru.h" // tune structures
#include "reorganize.h" // reorganize the geometry, cell, and wfc files
#include "directional.h" // estimate the directional H-bonds in liquid water
#include "dielectric.h" // compute the dielectric constants
#include "xsf.h" // compute the xsf files (for wave functions)
#include "density2D.h" // density related properties 
#include "HBs_near_PT.h" // calculate average accepted HBs of Hydroxide near PT. renxi 20200417
#include "qe_input.h" // transform every snapshots of a geo file into each qe input file. renxi 20200521
#include "PT_snapshot.h" // transform snapshots around PT into qe input and json files. renxi 20200615
#include "PT_snapshot2.h"
#include "fp_check.h"
#include "band_gap.h"
#include "incremental_pdf.h"
#include "special_msd.h"
#include "stress_average.h"
#include "mass_center.h"
#include "HB_angle.h"
#include "OH_movie.h"
#include "movie2.h"
#include "HB_stat.h"
#include "HB_stat2.h"
#include "HB_stat3.h"
#include "HB_stat4.h"
#include "Wan_centers_stat.h"
#include "dist2.h"
#include "incremental_pdf2.h"
#include "incremental_pdf3.h"
#include "HBs_near_AngularJump.h"
#include "HBs_near_DoubleDonor.h"
#include "HBs_near_TwistThird.h"
#include "HB_correlation.h"
#include "HB_correlation2.h"
#include "nonHB_correlation.h"
#include "nonHB_correlation2.h"
#include "nonHB_correlation3.h"
#include "HB_break.h"
#include "oho_angle.h"
#include "orientation_tcf.h"
#include "bdf_rcut1.h"
#include "first_shell_angle.h"
#include "Cation_free_energy.h"  // Jianchuan Liu add 2022-09-07
#include "OcationOAngle.h"  // Jianchuan Liu add 2022-09-07
#include "pos2pdb.h"  // Jianchuan Liu add 2022-09-07
int main(int argc, char **argv)
{
	//cout << "run D310" << endl;

#ifdef __MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &NPROC);
	MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
#endif

	//cout << "after initialization of MPI" << endl;

	stringstream ss;
	ss << "running" << RANK << ".log";
	ofs_running.open(ss.str().c_str());


	// read in the parameters.
	// INPUT has been generated in input.cpp
	INPUT.Init(argc,argv);
			
	//--------------------------------------------------
	// The program has several subroutines.
	// Each subroutine has its unit function.
	// I divide the functions into several categories.
	// Here are the explanations:
	//--------------------------------------------------
	// Category A: Structure related programs.
	// ext: Extend the cell periodically.
	// vacuum: Add vacuum to the input geometry and 
	//         output a new geometry.
	// pp: Plot the pseudopotentials.
	//--------------------------------------------------
	// Category B: Static structure properties.
	// pdf: Calcualte the pair distribution funnctions.
	// pdf2d: Calcualte the pair distribution funnctions on 2D plane.
	// pdf5: Calculate the pari distribution functions for the first few neighobours.
	// ssf: Calculate the static structure factors 
	//      in REAL SPACE.
	// vel: Calculate the distribution of velocities.
	// eprofile: Plot the profile of electron density.
	// iprofile: Plot the profile of ion density.
	//--------------------------------------------------
	// Category C: Dynamics structure properties.
	// dsf: dynamics structure factor.
	//--------------------------------------------------

	// Category A:
	if(INPUT.calculation == "ext"){Extend ext; ext.Routine();}
	else if(INPUT.calculation == "write"){Write write; write.Routine();}
	else if(INPUT.calculation == "tune_stru"){Tune_Stru ts; ts.Routine();}
	else if(INPUT.calculation == "vacuum"){Vacuum vacuum; vacuum.Routine();}
	else if(INPUT.calculation == "void"){Void voids; voids.Routine();}
	else if(INPUT.calculation == "pp"){Pseudo pp; pp.Routine();}
	else if(INPUT.calculation == "average_ions"){Average ave; ave.Routine();}
	else if(INPUT.calculation == "reorganize"){Reorganize RO; RO.Routine();}
	// Category B:
	else if(INPUT.calculation == "pdf"){PDF pdf; pdf.Routine();}
	else if(INPUT.calculation == "pdf2d"){PDF2d pdf2d; pdf2d.Routine();}
	else if(INPUT.calculation == "pdf5"){PDF5 pdf5; pdf5.Routine();}
	else if(INPUT.calculation == "ssf"){SSF ssf; ssf.Routine();}
	else if(INPUT.calculation == "ssf_selected"){SSF_Selected ssf2; ssf2.Routine();}
	else if(INPUT.calculation == "vel"){Vel vel; vel.Routine();}
	else if(INPUT.calculation == "eprofile"){Data3D data; data.Routine();}
	else if(INPUT.calculation == "iprofile"){Iprofile iprof; iprof.Routine();}
	else if(INPUT.calculation == "bdf"){BDF bdf; bdf.Routine();}
	else if(INPUT.calculation == "bdf_rcut"){BDF_Rcut bdf; bdf.Routine();}
	else if(INPUT.calculation == "ili"){ILI ili; ili.Routine();}
	else if(INPUT.calculation == "ili3d"){ILI_3D ili3; ili3.Routine();}
	else if(INPUT.calculation == "mdp"){MDP mdp; mdp.Routine();}
	else if(INPUT.calculation == "mdp2"){MDP2 mdp2; mdp2.Routine();}
	else if(INPUT.calculation == "mdp3"){MDP3 mdp3; mdp3.Routine();}
	else if(INPUT.calculation == "hbs"){HBs hbs; hbs.Routine();}
	else if(INPUT.calculation == "mj"){MJ mj; mj.Routine();}
	else if(INPUT.calculation == "top"){TOP top; top.Routine();}
	else if(INPUT.calculation == "pre"){PRE pre; pre.Routine();}
	else if(INPUT.calculation == "waterwire"){Waterwire ww; ww.Routine();}
	else if(INPUT.calculation == "waterwire2"){Waterwire2 ww2; ww2.Routine();}
	else if(INPUT.calculation == "ww_compress"){WW_Compress WWC; WWC.Routine();}
	else if(INPUT.calculation == "hyper"){Hyper hyper; hyper.Routine();}
	else if(INPUT.calculation == "msd"){MSD msd; msd.Routine();} 
	else if(INPUT.calculation == "msd_multiple"){MSD_Multiple msdm; msdm.Routine();} 
	else if(INPUT.calculation == "honeycutt"){Honeycutt ha; ha.Routine();} 
	// Category B2: 2D information
	else if(INPUT.calculation == "trajadj"){Trajadj tj; tj.Routine();}
	else if(INPUT.calculation == "xy_profile"){XY_Profile xy; xy.Routine();}
	else if(INPUT.calculation == "directional"){Directional direct; direct.Routine();}
	// Category C: Dynamic Properties
	else if(INPUT.calculation == "velcor"){VelCor velcor; velcor.Routine();}
	else if(INPUT.calculation == "powerspectra"){PowerSpectra ps; ps.Routine();}
	else if(INPUT.calculation == "isf"){ISF isf; isf.Routine();}
	else if(INPUT.calculation == "isf2"){ISF2 isf2; isf2.Routine();}//qianrui add 
	else if(INPUT.calculation == "dsf"){DSF dsf; dsf.Routine();}
	// Category D:
	else if(INPUT.calculation == "insert"){Insert inst; inst.Routine();}
	else if(INPUT.calculation == "movie"){Movie mov; mov.Routine();}
	else if(INPUT.calculation == "wannier"){Wannier wan; wan.Routine();} // compute dipole moments
	else if(INPUT.calculation == "wannier1"){Wannier1 wan1; wan1.Routine();} // analyze Wannier functions 
	else if(INPUT.calculation == "infrared"){Wannier wan; wan.infrared();} // compute infrared spectra
	else if(INPUT.calculation == "dielectric"){Dielectric die; die.Routine();} // compute infrared spectra
	else if(INPUT.calculation == "eig"){Eig eig; eig.Routine();}
	else if(INPUT.calculation == "dist"){Dist dist; dist.Routine();}
	else if(INPUT.calculation == "ele_conductivity"){Ele_Conductivity ele_con; ele_con.Routine();}//
	else if(INPUT.calculation == "elecond_contribute"){Elecond_contribute ele_contri; ele_contri.Routine();}
	else if(INPUT.calculation == "xsf"){XSF xsf; xsf.Routine();} // compute XSF related properties
	else if(INPUT.calculation == "density2d"){Density2D den; den.Routine();}
	else if(INPUT.calculation == "HBs_near_PT"){HBs_near_PT hbsn; hbsn.Routine();} // renxi added 20200417
	else if(INPUT.calculation == "qe_input"){qe_input qe; qe.Routine();}
	else if(INPUT.calculation == "PT_snapshot"){PT_snapshot ptss; ptss.Routine();}
	else if(INPUT.calculation == "PT_snapshot2"){PT_snapshot2 ptss2; ptss2.Routine();}
	else if(INPUT.calculation == "fp_check"){fp_check fpck; fpck.Routine();}
	else if(INPUT.calculation == "band_gap"){band_gap bg; bg.Routine();}
	else if(INPUT.calculation == "incremental_pdf"){incrementalPDF ip; ip.Routine();}
	else if(INPUT.calculation == "special_msd"){special_MSD sm; sm.Routine();}
	else if(INPUT.calculation == "stress_average"){stress_average sa; sa.Routine();}
        else if(INPUT.calculation == "mass_center"){mass_center mc; mc.Routine();}
	else if(INPUT.calculation == "hb_angle"){HB_angle hb; hb.Routine();}
	else if(INPUT.calculation == "oh_movie"){OH_movie hb; hb.Routine();}
	else if(INPUT.calculation == "movie2"){movie2 hb; hb.Routine();}
	else if(INPUT.calculation == "hb_stat"){HB_stat hb; hb.Routine();}
	else if(INPUT.calculation == "hb_stat2"){HB_stat2 hb; hb.Routine();}
	else if(INPUT.calculation == "hb_stat3"){HB_stat3 hb; hb.Routine();}
	else if(INPUT.calculation == "hb_stat4"){HB_stat4 hb; hb.Routine();}
	else if(INPUT.calculation == "hbs_near_angularjump"){HBs_near_AngularJump hb; hb.Routine();}
	else if(INPUT.calculation == "hbs_near_doubledonor"){HBs_near_DoubleDonor hb; hb.Routine();}
	else if(INPUT.calculation == "wan_centers_stat"){Wan_centers_stat hb; hb.Routine();}
	else if(INPUT.calculation == "dist2"){Dist2 hb; hb.Routine();}
	else if(INPUT.calculation == "incremental_pdf2"){incrementalPDF2 hb; hb.Routine();}
	else if (INPUT.calculation == "hbs_near_twistthird"){HBs_near_TwistThird tn; tn.Routine();}
	else if (INPUT.calculation == "hb_correlation"){HB_correlation tn; tn.Routine();}
	else if (INPUT.calculation == "hb_correlation2"){HB_correlation2 hb; hb.Routine();}
	else if (INPUT.calculation == "nonhb_correlation") {nonHB_correlation hb; hb.Routine();}
	else if (INPUT.calculation == "nonhb_correlation2") {nonHB_correlation2 hb; hb.Routine();}
	else if (INPUT.calculation == "nonhb_correlation3") {nonHB_correlation3 hb; hb.Routine();}
	else if (INPUT.calculation == "hb_break") {HB_break hb; hb.Routine();}
	else if (INPUT.calculation == "oho_angle") {oho_angle oho_angle; oho_angle.Routine();}
	else if(INPUT.calculation == "incremental_pdf3"){incrementalPDF3 hb; hb.Routine();}
	else if(INPUT.calculation == "orientation_tcf"){Orientation_TCF hb; hb.Routine();}
	else if(INPUT.calculation == "bdf_rcut1"){BDF_rcut1 hb; hb.Routine();}
	else if(INPUT.calculation == "first_shell_angle"){First_shell_angle hb; hb.Routine();}
    // free energy  CationFreeEnergy
	else if(INPUT.calculation == "CFE") {CationFreeEnergy cfe; cfe.Routine();}   // Jianchuan Liu add 2022-09-07
	// convert geometry to PDB format file. 
    else if(INPUT.calculation == "pos2pdb"){Pos2pdb pos2pdb; pos2pdb.Routine();} // Jianchuan Liu add 2022-09-07
	// O-Cation-O angle
	else if(INPUT.calculation == "o_c_o"){OcationOAngle oco; oco.Routine();}

	else
	{
		ofs_running << " calculation=" << INPUT.calculation << endl;
		QUIT("No 'calculation' available");
	}

	ofs_running << " --------------- " << endl;
	ofs_running << "      Finish     " << endl;
	ofs_running << " --------------- " << endl;

#ifdef __MPI
    MPI_Finalize();
#endif

	time_t time_end = std::time(nullptr);
	ofs_running << " HAVE A GREAT DAY! " << ctime(&time_end) << endl;

	return 0;
}
