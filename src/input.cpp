#include "input.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <string.h>

Input INPUT;


Input::Input() 
{
	bdf_nadj = 2;
	bdf_dtheta = 0.5;
	bdf_movie = 0;

	struf_dgx = 0.0;
	struf_dgy = 0.0;
	struf_dgz = 0.0;
	struf_ng = 0;	

	ext_1 = 0;
	ext_2 = 0;
	ext_3 = 0;

    geo_1 = 0;
	geo_2 = 0;
	geo_target = 0;

	nx = 0;
	ny = 0;
	nz = 0;
	x0 = -1.5;
	y0 = 4.5;
    //z0 = 0.0;
	dx = 0.05;
	dy = 0.05;
	//dz = 0.05; 
	dis_oc = -1;   // distace between cation and OH   // Jianchuan Liu add 2022-09-07
	model_oh = -1; // 0: Calculate the OH around Cation, 1: Calculate the OH not around Cation  // Jianchuan Liu add 2022-09-07
	// calculate the free energy  // Jianchuan Liu add 2022-09-07
	cfe_model = 1; // 1: X: distance of OH1-Cation, Y: distance of OH2-Cation
	               // 2: X: distance of OH1-Cation or OH2-Cation, Y: distance of OH-OH
				   // 3: X: avg distance of OH1-Cation and OH2-Cation, Y:  distance of OH-OH

	// for reorientation TCF Jianchuan Liu add 2022-09-07
	tcf_t = 6;  // ps 
	tcf_n = 1; // number of tcf needed
	tcf_t0 = 0.0; // starting point of tcf
	tcf_dt0 = 3.0; // difference between starting points
	tcf_dt = -0.0005; // unit is ps, time difference between two snapshots
    tcf_natom = -1;
    skip_frame = 1;


	triclinic = 0; // default is 0 for LAMMPS geometry file, mohan add 2015-05-15	
	hindex = 0;
	ele1 = "none";
	ele2 = "none";
	bdf_center = "none";
	z1 = -10000;
	z0 = -10000;
    rho_ref = 0.0;
	delta = 0.0;
    ion_analysis = false;
	rcut1 = 10; //Angstroms
    id1 = "none";
    id2 = "none";
	id3 = "none";
	ntry = 50;
	dz = 0;
	ili_file = "none";
	func = 1;
	func_b = 1;
	func_c = 1;
	func_d = 1;
	func_e = 1;
	factor = 1.0;
	within = 0.01;
	maxiter = 100;
	dr = 0.01;
	rcut = 6.00;
	ref_rho = 0.016;
    //upper_z = 25.5; // for hexane+water
    upper_z = 10000; // mohan update 2017-09-19 
	u1 = 20;
	u2 = 20;
	u3 = 20; // mohan add 2016-11-26
	a0 = -2.0;
	a1 = 0.01;
    mdp0 = -2.0;
	only_hydroxide=false;
	rcut_oo = 3.50; // unit is Angstroms
	rcut_oh = 1.24; // unit is Angstroms
	rcut_ch = 1.24; // unit is Angstroms
	rcut_clo = 3.70; // unit is Angstroms
	acut_hoo = 30; // cutoff of angle in degrees
	d_angle = 0.5; // delta angle in degrees
	geo_format = 1;
	system = "none";

	snatom = -1;
	satom = new int[1];

	nacc = -1;
    ndon = -1;
	// for wannier
	nbands = 0;
	wannier_file = "none";
    eig_file = "none";
	pos_ili_file = "none";
	shift = 0.0;
	cell_file = "none";

	// for trajadj
	ntzone = 0;

	// for iprof
	iprof_nr = 0;
	iprof_b = 0.0;
    nbin = 100; // default value

	celldm1 = 0.0; celldm2 = 0.0; celldm3 = 0.0;
	e11 = 0.0; e12 = 0.0; e13 = 0.0;
	e21 = 0.0; e22 = 0.0; e23 = 0.0;
	e31 = 0.0; e32 = 0.0; e33 = 0.0;

	length_unit = "bohr";

	// HA_index
	HA_pdf = false;
	HA_nsn = 5;
	HA_nsb = 5;

	movement_x=0.0;
	movement_y=0.0;
	movement_z=0.0;

	npoints = 1;
	dt_snapshots = -1; // unit is ps 
	dt_max = 0.5;
    ps_dt = 0.005; // unit is ps

	msd_n = 1; // number of MSD needed
	msd_t0 = 0.0; // starting point of MSD
	msd_t = 50.0; // unit is ps
	msd_dt0 = 3.0; // difference between starting points
	msd_dt = -0.0005; // unit is ps, time difference between two snapshots, mohan updated 2019-04-05
	msd_natom=1; //number of atoms
	msd_stokes=false; // the hydrodynamics diffusion of water ions

	// for infrared spectra
	tcor = 1000;

	mu = 0.0;

	temperature = 330; // unit is K

	// number of neighbours.
	neighbours = 5;

	geo_out = "result.dat";

	cartesian = true;

	force_file = "none";
	
	outdir = "none";

	// for PT_snapshot renxi 20200615
	nPT = 0;
	geo_target_aft = 0;

	// for PT_snapshot2 renxi 20200628
	output_interval = 0;
	nsnapshot = 0;

	// for band_gap renxi 20200912

	HO = 0;
	LU = 0;

	// for incremental pdf renxi 20200919
	nshell = -1;

	// for special_MSD renxi 20200924
	ia_select = -1;

	// for stress average renxi 20200926
	stress_geometry = "none";

	// for water_within_distance renxi 20210707

	radius_inner = -1;
	radius_outer = -1;

	vmax = -1;
	dv = -1;
	stay_tmax = -1;
	stay_dt = -1;

	vel_file = "none";
	dtheta = -1;

	nHB_max = 8;

	HB_defination = 1;
	HB_defination_k = 1000;
	HB_defination_b = -1;

	relax_time = 0.5;
	dq = 0.01;
	theta_min = 0;
	r_max = 100;
	r_min = 0;

	pdf_nstd = -1;
	abacus_version = "old";

	n_recorded_water = 0;
	Oindex = new int[10];

	non_return = 0;
	upper_time = 0;
	lower_time = 0;
}

Input::~Input() 
{
	delete[] satom;
}

void Input::Init(const string &fn, const int &argc)
{
	time_t time_start = std::time(NULL);
	ofs_running << " ------------------------------------------------------------------------------------" << endl;
	ofs_running << "              WELCOME TO D310 " << ctime(&time_start) << endl;
	ofs_running << " THIS IS A PROGRAM FOR PERFORMING ANALYSIS FOR MOLECULAR DYNAMICS" << endl;
	ofs_running << " AUTHOR: MOHAN CHEN, LAST UPDATE: 2017-07" << endl;

	ofs_running << setiosflags(ios::right);
	ofs_running << setiosflags(ios::left);
	ofs_running << setiosflags(ios::left);

	ofs_running << " READING from INPUT file" << endl;


    this->Default();

    this->Read(fn);

    this->Check();

    time_t  time_now = time(NULL);

	return;
}

void Input::Default(void)
{
	step_interval_dynamics = -1;

	geo_interval = -1;
    geo_ignore = -1;

	format3D = 1;

	isf_dg = -1.0;
	isf_ng = -1;
	isf_nt1=0;//qianrui
    isf_nt2=0;//qianrui
    isf_m1=0;//qianrui
    isf_m2=0;//qianrui
    isf_m3=0;//qianrui

	natom1 = 0;
	natom2 = 0;
	natom3 = 0;
	natom4 = 0;

	id1 = "empty";
	id2 = "empty";
	id3 = "empty";
	id4 = "empty";
	theta = 90;
    return;
}

void Input::Read(const string &fn)
{
    ifstream ifs(fn.c_str(), ios::in);	// "in_datas/input_parameters"

    if (!ifs) 
	{
		cout << " Can't find the INPUT file." << endl;
		exit(0);
	}

    ifs.clear();
    ifs.seekg(0);

    char word[80];
    char word1[80];
    int ierr = 0;

    //ifs >> setiosflags(ios::uppercase);
    ifs.rdstate();

    while (ifs.good())
    {
        ifs >> word1;
        strtolower(word1, word);
//		cout << " word = " << word << endl;

//----------------------------------------------------------
// main parameters
//----------------------------------------------------------
		// >>> General <<<
        if (strcmp("calculation", word) == 0) read_value(ifs, calculation);
        else if (strcmp("ext_1", word) == 0) read_value(ifs, ext_1);
        else if (strcmp("ext_2", word) == 0) read_value(ifs, ext_2);
        else if (strcmp("ext_3", word) == 0) read_value(ifs, ext_3);
        else if (strcmp("geo_in", word) == 0) read_value(ifs, geo_in);
        else if (strcmp("geo_in_type", word) == 0) read_value(ifs, geo_in_type);
        else if (strcmp("geo_out", word) == 0) read_value(ifs, geo_out);
        else if (strcmp("geo_out_type", word) == 0) read_value(ifs, geo_out_type);
	else if (strcmp("geo_directory", word) == 0) read_value(ifs, geo_directory);
	else if (strcmp("geo_1", word) == 0) read_value(ifs, geo_1);
	else if (strcmp("geo_2", word) == 0) read_value(ifs, geo_2);
	else if (strcmp("geo_target", word) == 0) read_value(ifs, geo_target);
	else if (strcmp("geo_interval", word) == 0) read_value(ifs, geo_interval);
	else if (strcmp("geo_ignore", word) == 0) read_value(ifs, geo_ignore);
	else if (strcmp("geo_format", word) == 0) read_value(ifs, geo_format);
	else if (strcmp("ntype", word) == 0) read_value(ifs, ntype);
        else if (strcmp("natom", word) == 0) read_value(ifs, natom);
        else if (strcmp("cartesian", word) == 0) read_value(ifs, cartesian); //mohan add 2014-04-04
		// >>> Function 2 <<<
        else if (strcmp("dr", word) == 0) read_value(ifs, dr);
        else if (strcmp("rcut", word) == 0) read_value(ifs, rcut);
        else if (strcmp("pdf_z0", word) == 0) read_value(ifs, pdf_z0);
        else if (strcmp("pdf_z1", word) == 0) read_value(ifs, pdf_z1);
        else if (strcmp("rho_ref", word) == 0) read_value(ifs, rho_ref);
        else if (strcmp("delta", word) == 0) read_value(ifs, delta);
        else if (strcmp("struf_dgx", word) == 0) read_value(ifs, struf_dgx);
        else if (strcmp("struf_dgy", word) == 0) read_value(ifs, struf_dgy);
        else if (strcmp("struf_dgz", word) == 0) read_value(ifs, struf_dgz);
        else if (strcmp("struf_ng", word) == 0) read_value(ifs, struf_ng);
        else if (strcmp("ssf_out", word) == 0) read_value(ifs, ssf_out);
		// >>> Function 3 <<<
        else if (strcmp("vel_in", word) == 0) read_value(ifs, vel_in);
        else if (strcmp("vel_out", word) == 0) read_value(ifs, vel_out);
        else if (strcmp("ndv", word) == 0) read_value(ifs, ndv);
		// >>> Function 4 <<<
        else if (strcmp("vacuum_x1", word) == 0) read_value(ifs, vacuum_x1);
        else if (strcmp("vacuum_x2", word) == 0) read_value(ifs, vacuum_x2);
        else if (strcmp("vacuum_y1", word) == 0) read_value(ifs, vacuum_y1);
        else if (strcmp("vacuum_y2", word) == 0) read_value(ifs, vacuum_y2);
        else if (strcmp("vacuum_z1", word) == 0) read_value(ifs, vacuum_z1);
        else if (strcmp("vacuum_z2", word) == 0) read_value(ifs, vacuum_z2);
		// >>> Function 5 <<<
        else if (strcmp("direction", word) == 0) read_value(ifs, direction);
        else if (strcmp("data_in", word) == 0) read_value(ifs, data_in);
        else if (strcmp("data_out", word) == 0) read_value(ifs, data_out);
		// >>> Function : pseudopotential <<<
        else if (strcmp("pseudo_z", word) == 0) read_value(ifs, pseudo_z);
        else if (strcmp("pseudo_type", word) == 0) read_value(ifs, pseudo_type);
        else if (strcmp("pseudo_in", word) == 0) read_value(ifs, pseudo_in);
        else if (strcmp("pseudo_out", word) == 0) read_value(ifs, pseudo_out);
		// >>> Function : iprofile
        else if (strcmp("iprof_nr", word) == 0) read_value(ifs, iprof_nr);
        else if (strcmp("iprof_b", word) == 0) read_value(ifs, iprof_b);
        else if (strcmp("nbin", word) == 0) read_value(ifs, nbin);
		// >>> Function : eprofile
        else if (strcmp("format3d", word) == 0) read_value(ifs, format3D);
		// >>> Function : dynamics structure factor
        else if (strcmp("dsf_out", word) == 0) read_value(ifs, dsf_out);
        else if (strcmp("dsf_dt", word) == 0) read_value(ifs, dsf_dt);
        else if (strcmp("dsf_neqi", word) == 0) read_value(ifs, dsf_neqi);
		// >>> Function : velocity correlation function
        else if (strcmp("velcor_in_type", word) == 0) read_value(ifs, velcor_in_type);
        else if (strcmp("velcor_directory", word) == 0) read_value(ifs, velcor_directory);
		else if (strcmp("velcor_1", word) == 0) read_value(ifs, velcor_1);
		else if (strcmp("velcor_2", word) == 0) read_value(ifs, velcor_2);
		else if (strcmp("velcor_neqi", word) == 0) read_value(ifs, velcor_neqi);
		else if (strcmp("velcor_out", word) == 0) read_value(ifs, velcor_out);
		else if (strcmp("velcor_atom", word) == 0) read_value(ifs, velcor_atom);
		else if (strcmp("step_interval_dynamics", word) == 0) read_value(ifs, step_interval_dynamics);
		// >>> Function : power spectra
		else if (strcmp("ps_nv", word) == 0) read_value(ifs, ps_nv);
		else if (strcmp("ps_dw", word) == 0) read_value(ifs, ps_dw);
		else if (strcmp("ps_nw", word) == 0) read_value(ifs, ps_nw);
		else if (strcmp("ps_in", word) == 0) read_value(ifs, ps_in);
		else if (strcmp("ps_out", word) == 0) read_value(ifs, ps_out);
		else if (strcmp("ps_dt", word) == 0) read_value(ifs, ps_dt);
		// >>> Function : intermediate scattering function
		else if (strcmp("isf_target_q", word) == 0) read_value(ifs, isf_target_q);
		else if (strcmp("isf_nconfig", word) == 0) read_value(ifs, isf_nconfig);
		else if (strcmp("isf_ncorrelation", word) == 0) read_value(ifs, isf_ncorrelation);
		else if (strcmp("isf_dcorrelation", word) == 0) read_value(ifs, isf_dcorrelation);
		else if (strcmp("isf_outfile", word) == 0) read_value(ifs, isf_outfile);
		else if (strcmp("natom1", word) == 0) read_value(ifs, natom1);
		else if (strcmp("natom2", word) == 0) read_value(ifs, natom2);
		else if (strcmp("natom3", word) == 0) read_value(ifs, natom3);
		else if (strcmp("natom4", word) == 0) read_value(ifs, natom4);
		else if (strcmp("isf_dg", word) == 0) read_value(ifs, isf_dg);
		else if (strcmp("isf_ng", word) == 0) read_value(ifs, isf_ng);
		else if (strcmp("isf_nt1", word) == 0) read_value(ifs, isf_nt1);//qianrui
        else if (strcmp("isf_nt2", word) == 0) read_value(ifs, isf_nt2);//qianrui
        else if (strcmp("isf_m1", word) == 0) read_value(ifs, isf_m1);//qianrui
        else if (strcmp("isf_m2", word) == 0) read_value(ifs, isf_m2);//qianrui
        else if (strcmp("isf_m3", word) == 0) read_value(ifs, isf_m3);//qianrui
		else if (strcmp("isf_ngx", word) == 0) read_value(ifs, isf_ngx);//qianrui
        else if (strcmp("isf_ngy", word) == 0) read_value(ifs, isf_ngy);//qianrui
        else if (strcmp("isf_ngz", word) == 0) read_value(ifs, isf_ngz);//qianrui
        else if (strcmp("isf_dgx", word) == 0) read_value(ifs, isf_dgx);//qianrui
        else if (strcmp("isf_dgy", word) == 0) read_value(ifs, isf_dgy);//qianrui
        else if (strcmp("isf_dgz", word) == 0) read_value(ifs, isf_dgz);//qianrui

		else if (strcmp("isf_config_start", word) == 0) read_value(ifs, isf_config_start);
		// >>> Function : bond-angle distribution functions
		else if (strcmp("bdf_nadj", word) == 0) read_value(ifs, bdf_nadj);
		else if (strcmp("bdf_dtheta", word) == 0) read_value(ifs, bdf_dtheta);
		else if (strcmp("bdf_movie", word) == 0) read_value(ifs, bdf_movie);
		else if (strcmp("bdf_rcut", word) == 0) read_value(ifs, bdf_rcut);
		else if (strcmp("bdf_center", word) == 0) read_value(ifs, bdf_center);
		else if (strcmp("rcut1", word) == 0) read_value(ifs, rcut1);
		else if (strcmp("z0", word) == 0) read_value(ifs, z0);
		else if (strcmp("z1", word) == 0) read_value(ifs, z1);
		// >>> Function : movement z
		else if (strcmp("movement_x", word) == 0) read_value(ifs, movement_x);
		else if (strcmp("movement_y", word) == 0) read_value(ifs, movement_y);
		else if (strcmp("movement_z", word) == 0) read_value(ifs, movement_z);
		// >>> Function : insert atoms
		else if (strcmp("natom_new", word) == 0) read_value(ifs, natom_new);
		else if (strcmp("element_new", word) == 0) read_value(ifs, element_new);
		else if (strcmp("min_dis", word) == 0) read_value(ifs, min_dis);
        // >>> Function : triclinic
		else if (strcmp("triclinic", word) == 0) read_value(ifs, triclinic); //mohan add 2015-06-16
		// >>> Function : QE geometry
		else if (strcmp("cell_file", word) == 0) read_value(ifs, cell_file); //mohan add 2016-12-22
		else if (strcmp("celldm1", word) == 0) read_value(ifs, celldm1); //mohan add 2016-08-02
		else if (strcmp("celldm2", word) == 0) read_value(ifs, celldm2); //mohan add 2016-08-02
		else if (strcmp("celldm3", word) == 0) read_value(ifs, celldm3); //mohan add 2016-08-02
		else if (strcmp("e11", word) == 0) read_value(ifs, e11); //mohan add 2018-02-23
		else if (strcmp("e12", word) == 0) read_value(ifs, e12); //mohan add 2018-02-23
		else if (strcmp("e13", word) == 0) read_value(ifs, e13); //mohan add 2018-02-23
		else if (strcmp("e21", word) == 0) read_value(ifs, e21); //mohan add 2018-02-23
		else if (strcmp("e22", word) == 0) read_value(ifs, e22); //mohan add 2018-02-23
		else if (strcmp("e23", word) == 0) read_value(ifs, e23); //mohan add 2018-02-23
		else if (strcmp("e31", word) == 0) read_value(ifs, e31); //mohan add 2018-02-23
		else if (strcmp("e32", word) == 0) read_value(ifs, e32); //mohan add 2018-02-23
		else if (strcmp("e33", word) == 0) read_value(ifs, e33); //mohan add 2018-02-23
		else if (strcmp("length_unit", word) == 0) read_value(ifs, length_unit); //mohan add 2016-08-02
		else if (strcmp("id1", word) == 0) read_value(ifs, id1); //mohan add 2016-08-02
		else if (strcmp("id2", word) == 0) read_value(ifs, id2); //mohan add 2016-08-02
		else if (strcmp("id3", word) == 0) read_value(ifs, id3); //mohan add 2016-08-02
		else if (strcmp("id4", word) == 0) read_value(ifs, id4); //mohan add 2016-08-02
		else if (strcmp("ele1", word) == 0) read_value(ifs, ele1); //mohan add 2016-08-02
		else if (strcmp("ele2", word) == 0) read_value(ifs, ele2); //mohan add 2016-08-02
		else if (strcmp("ion_analysis", word) == 0) read_value(ifs, ion_analysis); //mohan add 2016-10-04
		// for iprof
		else if (strcmp("hindex", word) == 0) read_value(ifs, hindex); //mohan add 2015-06-16
		// for ili
		else if (strcmp("nx", word) == 0) read_value(ifs, nx); //mohan add 2016-10-11 
		else if (strcmp("ny", word) == 0) read_value(ifs, ny); //mohan add 2016-10-11 
		else if (strcmp("nz", word) == 0) read_value(ifs, nz); //mohan add 2016-10-11 
		else if (strcmp("d", word) == 0) read_value(ifs, d); //mohan add 2016-10-11 
		else if (strcmp("zeta", word) == 0) read_value(ifs, zeta); //mohan add 2016-10-11 
		else if (strcmp("ntry", word) == 0) read_value(ifs, ntry); //mohan add 2016-10-11 
		else if (strcmp("dz", word) == 0) read_value(ifs, dz); //mohan add 2016-10-11 
		else if (strcmp("ref_rho", word) == 0) read_value(ifs, ref_rho); //mohan add 2016-10-11 
		else if (strcmp("upper_z", word) == 0) read_value(ifs, upper_z); //mohan add 2016-10-11 
		else if (strcmp("dx", word) == 0) read_value(ifs, dx); //mohan add 2017-01-03 
		else if (strcmp("dy", word) == 0) read_value(ifs, dy); //mohan add 2017-01-03 
//		else if (strcmp("dz", word) == 0) read_value(ifs, dz); //mohan add 2018-08-20 
		else if (strcmp("x0", word) == 0) read_value(ifs, x0); //mohan add 2017-01-03 
		else if (strcmp("y0", word) == 0) read_value(ifs, y0); //mohan add 2017-01-03 
//		else if (strcmp("z0", word) == 0) read_value(ifs, z0); //mohan add 2018-08-20 
		// for mdp that based on ili
		else if (strcmp("ili_file", word) == 0) read_value(ifs, ili_file); //mohan add 2016-10-11 
		else if (strcmp("u1", word) == 0) read_value(ifs, u1); //mohan add 2016-10-11 
		else if (strcmp("u2", word) == 0) read_value(ifs, u2); //mohan add 2016-10-11 
		else if (strcmp("u3", word) == 0) read_value(ifs, u3); //mohan add 2016-10-11 
		else if (strcmp("only_hydroxide", word) == 0) read_value(ifs, only_hydroxide); //mohan add 2016-10-11 
		else if (strcmp("a0", word) == 0) read_value(ifs, a0); //mohan add 2016-10-11 
		else if (strcmp("a1", word) == 0) read_value(ifs, a1); //mohan add 2016-10-11 
		else if (strcmp("mdp0", word) == 0) read_value(ifs, mdp0); //mohan add 2016-10-11 
		// for water and hydrogen bonds
		else if (strcmp("rcut_oo", word) == 0) read_value(ifs, rcut_oo); //mohan add 2016-10-21
		else if (strcmp("rcut_oh", word) == 0) read_value(ifs, rcut_oh); //mohan add 2016-10-21
		else if (strcmp("rcut_ch", word) == 0) read_value(ifs, rcut_ch); //mohan add 2019-03-17
		else if (strcmp("rcut_clo", word) == 0) read_value(ifs, rcut_clo); //mohan add 2016-10-21
		else if (strcmp("acut_hoo", word) == 0) read_value(ifs, acut_hoo); //mohan add 2016-10-21
		else if (strcmp("d_angle", word) == 0) read_value(ifs, d_angle); //mohan add 2016-10-21
		else if (strcmp("system", word) == 0) read_value(ifs, system); //mohan add 2016-10-21
		// HA
		else if (strcmp("ha_pdf", word) == 0) read_value(ifs, HA_pdf); //mohan add 2017-01-15
		else if (strcmp("ha_nsn", word) == 0) read_value(ifs, HA_nsn); //mohan add 2017-01-15
		else if (strcmp("ha_nsb", word) == 0) read_value(ifs, HA_nsb); //mohan add 2017-01-15
		// for general purposes
		else if (strcmp("factor", word) == 0) read_value(ifs, factor); //mohan add 2016-10-31
		else if (strcmp("within", word) == 0) read_value(ifs, within); //mohan add 2016-10-31
		else if (strcmp("maxiter", word) == 0) read_value(ifs, maxiter); //mohan add 2016-10-31
		else if (strcmp("func", word) == 0) read_value(ifs, func); //mohan add 2016-10-31
		else if (strcmp("func_b", word) == 0) read_value(ifs, func_b); //mohan add 2016-10-31
		else if (strcmp("func_c", word) == 0) read_value(ifs, func_c); //mohan add 2017-04-13
		else if (strcmp("func_d", word) == 0) read_value(ifs, func_d); //mohan add 2017-12-24
		else if (strcmp("func_e", word) == 0) read_value(ifs, func_e);
		else if (strcmp("npoints", word) == 0) read_value(ifs, npoints); //mohan add 2016-10-31
		else if (strcmp("dt_snapshots", word) == 0) read_value(ifs, dt_snapshots); //mohan add 2018-07-05
		// for infrared spectra
		else if (strcmp("tcor", word) == 0) read_value(ifs, tcor); // mohan add 2018-01-18 
		else if (strcmp("mu", word) == 0) read_value(ifs, mu); // mohan add 2018-01-18 
		// number of neighbours
		else if (strcmp("neighbours", word) == 0) read_value(ifs, neighbours); // mohan add 2018-04-09 
		// dipole file
        else if (strcmp("dipole_file", word) == 0) read_value(ifs, dipole_file);
        else if (strcmp("vdipole_file", word) == 0) read_value(ifs, vdipole_file);
        else if (strcmp("temperature", word) == 0) read_value(ifs, temperature); // mohan add 2019-02-07
		else if (strcmp("force_file", word) == 0) read_value(ifs, force_file);
		else if (strcmp("pdf_nstd", word) == 0) read_value(ifs, pdf_nstd);
		// cation free energy
		else if (strcmp("dis_oc", word) == 0) read_value(ifs, dis_oc); //Jianchuan Liu add 2021-11-05 
		else if (strcmp("cfe_model", word) == 0) read_value(ifs, cfe_model); //Jianchuan Liu add 2022-03-14 
		// TCF 
		else if (strcmp("tcf_t", word) == 0) read_value(ifs, tcf_t); //Jianchuan Liu add 2022-01-10 
		else if (strcmp("tcf_t0", word) == 0) read_value(ifs, tcf_t0); //Jianchuan Liu add 2022-01-10 
		else if (strcmp("tcf_dt0", word) == 0) read_value(ifs, tcf_dt0); //Jianchuan Liu add 2022-01-10 
		else if (strcmp("tcf_n", word) == 0) read_value(ifs, tcf_n); //Jianchuan Liu add 2022-01-10 
		else if (strcmp("tcf_dt", word) == 0) read_value(ifs, tcf_dt); //Jianchuan Liu add 2022-01-10
		else if (strcmp("tcf_natom", word) == 0) read_value(ifs, tcf_natom); //Jianchuan Liu add 2022-01-10
		else if (strcmp("skip_frame", word) == 0) read_value(ifs, skip_frame); //Jianchuan Liu add 2022-01-10
		else if (strcmp("model_oh", word) == 0) read_value(ifs, model_oh); //Jianchuan Liu add 2021-11-09
		else if (strcmp("ion_analysis", word) == 0) read_value(ifs, ion_analysis); //mohan add 2016-10-04
		// for movie
		else if (strcmp("snatom", word) == 0){ 
				ifs >> snatom;
				cout << "snatom=" << snatom << endl;
				assert(snatom>0);
				delete[] satom;
				this->satom = new int[snatom]();
				for(int i=0; i<snatom; ++i)
				{
					ifs >> satom[i];
					cout << "readin satom=" << satom[i] << endl;
				}
        		ifs.ignore(150, '\n');
				cout << "Read snatom done." << endl;
				// still has bugs! Read twice. mohan note 2016-10-06
		}
        else if (strcmp("nacc", word) == 0) read_value(ifs, nacc);
        else if (strcmp("ndon", word) == 0) read_value(ifs, ndon);
		// for wannier
        else if (strcmp("nbands", word) == 0) read_value(ifs, nbands);
        else if (strcmp("wannier_file", word) == 0) read_value(ifs, wannier_file);
        else if (strcmp("eig_file", word) == 0) read_value(ifs, eig_file);
        else if (strcmp("pos_ili_file", word) == 0) read_value(ifs, pos_ili_file);
        else if (strcmp("shift", word) == 0) read_value(ifs, shift);
 		// trajadj
		else if (strcmp("ntzone", word) == 0) read_value(ifs, ntzone);
		// for mean square displacements (MSD)
		else if (strcmp("msd_n", word) == 0) read_value(ifs, msd_n);
		else if (strcmp("msd_t0", word) == 0) read_value(ifs, msd_t0);
		else if (strcmp("msd_t", word) == 0) read_value(ifs, msd_t);
		else if (strcmp("msd_dt0", word) == 0) read_value(ifs, msd_dt0);
		else if (strcmp("msd_dt", word) == 0) read_value(ifs, msd_dt);
		else if (strcmp("msd_natom", word) == 0) read_value(ifs, msd_natom);
		else if (strcmp("msd_stokes", word) == 0) read_value(ifs, msd_stokes);
		else if (strcmp("outdir", word) == 0) read_value(ifs, outdir);
		else if (strcmp("npt", word) == 0) read_value(ifs, nPT);
		else if (strcmp("geo_target_aft", word) == 0) read_value(ifs, geo_target_aft);
		else if (strcmp("output_interval", word) == 0) read_value(ifs, output_interval);
		else if (strcmp("nsnapshot", word) == 0) read_value(ifs, nsnapshot);
		else if (strcmp("ho", word) == 0) read_value(ifs, HO);
		else if (strcmp("lu", word) == 0) read_value(ifs, LU);
		else if (strcmp("nshell", word) == 0) read_value(ifs, nshell);
		else if (strcmp("ia_select", word) == 0) read_value(ifs, ia_select);
		else if (strcmp("radius_inner", word) == 0) read_value(ifs, radius_inner);
		else if (strcmp("radius_outer", word) == 0) read_value(ifs, radius_outer);
		else if (strcmp("vmax", word) == 0) read_value(ifs, vmax);
		else if (strcmp("dv", word) == 0) read_value(ifs, dv);
		else if (strcmp("stay_tmax", word) == 0) read_value(ifs, stay_tmax);
		else if (strcmp("stay_dt", word) == 0) read_value(ifs, stay_dt);
		else if (strcmp("vel_file", word) == 0) read_value(ifs, vel_file);
		else if (strcmp("dtheta", word) == 0) read_value(ifs, dtheta);
		else if (strcmp("nhb_max", word) == 0) read_value(ifs, nHB_max);
		else if (strcmp("hb_defination", word) == 0) read_value(ifs, HB_defination);
		else if (strcmp("hb_defination_k", word) == 0) read_value(ifs, HB_defination_k);
		else if (strcmp("hb_defination_b", word) == 0) read_value(ifs, HB_defination_b);
		else if (strcmp("relax_time", word) == 0) read_value(ifs, relax_time);
		else if (strcmp("theta", word) == 0) read_value(ifs, theta);
		else if (strcmp("theta_min", word) == 0) read_value(ifs, theta_min);
		else if (strcmp("r_max", word) == 0) read_value(ifs, r_max);
		else if (strcmp("r_min", word) == 0) read_value(ifs, r_min);
		else if (strcmp("dq", word) == 0) read_value(ifs, dq);
		else if (strcmp("abacus_version", word) == 0) read_value(ifs, abacus_version);
		else if (strcmp("n_recorded_water", word) == 0) read_value(ifs, n_recorded_water);
		else if (strcmp("dt_max", word) == 0) read_value(ifs, dt_max);
		else if (strcmp("non_return", word) == 0) read_value(ifs, non_return);
		else if (strcmp("upper_time", word) == 0) read_value(ifs, upper_time);
		else if (strcmp("lower_time", word) == 0) read_value(ifs, lower_time);
		else if (strcmp("oindex", word) == 0)
		{
			for (int io=0; io<n_recorded_water-1; io++)
			{
				ifs >> Oindex[io];
			}
			read_value(ifs, Oindex[n_recorded_water-1]);
		}
		// >>> Function : new
        else
        {
            cout << " THE PARAMETER NAME '" << word
               << "' IS NOT USED!" << endl;
            ifs.ignore(150, '\n');
        }

        ifs.rdstate();
        if (ifs.eof() != 0)
        {
			break;
        }
        else if (ifs.bad() != 0)
        {
			cout << " Bad input parameters. " << endl;
			exit(0);
        }
        else if (ifs.fail() != 0)
        {
			cout << " word = " << word << endl;
			cout << " Fail to read parameters. " << endl; 
            ifs.clear();
			exit(0);
        }
        else if (ifs.good() == 0)
        {
			break;
        }
    }

    return;
}//end read_parameters


void Input::Check(void)
{
	if( vacuum_x1 < 0.0 ) QUIT("vacuum_x1 < 0.0"); 
	if( vacuum_x2 < 0.0 ) QUIT("vacuum_x2 < 0.0"); 
	if( vacuum_y1 < 0.0 ) QUIT("vacuum_y1 < 0.0"); 
	if( vacuum_y2 < 0.0 ) QUIT("vacuum_y2 < 0.0"); 
	if( vacuum_z1 < 0.0 ) QUIT("vacuum_z1 < 0.0"); 
	if( vacuum_z2 < 0.0 ) QUIT("vacuum_z2 < 0.0"); 

	ofs_running << " ------------------------------------------------------" << endl;
	if(calculation=="ssf") cout << " Category B: Static Structure Factor. " << endl;
	else if(calculation=="dsf") cout << " Category C: Dynamics Structure Factor. " << endl;
	else if(calculation=="velcor") cout << " Category C: Velocity Auto-Correlation Function. " << endl;
	ofs_running << " ------------------------------------------------------" << endl;

    return;
}


void Input::readbool(ifstream &ifs, bool &var)
{
    string str;
    ifs >> str;
    if (str == "true")
    {
        var = true;
    }
    else
    {
        var = false;
    }
    ifs.ignore(100, '\n');
    return;
}

void Input::strtolower(char *sa, char *sb)
{
    char c;
    int len = strlen(sa);
    for (int i = 0; i < len; i++)
    {
        c = sa[i];
        sb[i] = tolower(c);
    }
    sb[len] = '\0';
}

