#include "qe_input.h"
#include "input.h"

qe_input::qe_input(){

}

qe_input::~qe_input(){

}

void qe_input::Routine()
{
	ifstream ifs(INPUT.geo_directory);
	int snapshot_index;
	double snapshot_time;
	double x, y, z;
	for(int igeo = INPUT.geo_1; igeo <= INPUT.geo_2; igeo++)
	{
		ofstream ofs("./input_file/" + to_string(igeo) + "test.in");
		ofs << setprecision(14);
		ofs << "&CONTROL" << endl;
		ofs << "  title = \"water 63mol+H+ PBE 330K\"," << endl;
		ofs << "  calculation = 'cp-wf'," << endl;
		ofs << "  restart_mode = \"from_scratch\"" << endl;
		ofs << "  prefix = \"h+scan\"," << endl;
		ofs << "  pseudo_dir = '" << INPUT.pseudo_in << "'," << endl;
		ofs << "  outdir = '" << INPUT.outdir << "'," << endl;
		ofs << "  nstep = 1," << endl;
		ofs << "  isave = 1," << endl;
		ofs << "  dt = 2.D0," << endl;
		ofs << "  tstress = .true.," << endl;
		ofs << "  tprnfor = .true.," << endl;
		ofs << "  iprint = 1," << endl;
		ofs << "  ekin_conv_thr = 1.D-6," << endl;
		ofs << "  max_seconds = 43000," << endl;
		
		ofs << "/" << endl;
		ofs << "&SYSTEM" << endl;
		ofs << "  ibrav     = 1," << endl;
		ofs << "  celldm(1) = 23.5170," << endl;
		ofs << "  nat       = 190," << endl;
		ofs << "  ntyp      = 2," << endl;
		ofs << "  ecutwfc   = 85.00D0," << endl;
		ofs << "  ecfixed   = 130.D0, ! effective cutoff value" << endl;
		ofs << "  q2sigma   = 15.D0," << endl;
		ofs << "  qcutz     = 200.D0," << endl;
		ofs << "  input_dft = 'scan'," << endl;
		ofs << "  tot_charge    = +1," << endl;
		
		ofs << "/" << endl << "&ELECTRONS" << endl;
		ofs << "  electron_dynamics     = 'damp'," << endl;
		ofs << "  electron_damping      = 0.20," << endl;
		ofs << "  emass                 = 100.D0," << endl;
		ofs << "  emass_cutoff          = 25.D0," << endl;
		ofs << "  ortho_max             = 800," << endl;
		ofs << "  ortho_eps             = 1.D-7," << endl;
		ofs << "  electron_maxstep      = 500," << endl;

		ofs << "/" << endl << "IONS" << endl;
		ofs << "  ion_dynamics = \"none\"" << endl;
		ofs << "  ion_temperature = 'not_controlled'" << endl;
		ofs << "  fnosep           = 60," << endl;
		ofs << "  tempw            = 330," << endl;
		ofs << "  nhpcl            = 4," << endl;
		ofs << "  nhptyp	   = 2," << endl;
		ofs << "  ndega            = -3," << endl;
		ofs << "  ion_radius(1)    = 1.4," << endl;
		ofs << "  ion_radius(2)    = 1.4," << endl;

		ofs << "/" << endl << "&WANNIER" << endl;
		ofs << "  nit    =  60," << endl;
		ofs << "  calwf  =   3," << endl;
		ofs << "  tolw   =   1.D-6, nsteps =   50," << endl;
		ofs << "  adapt  = .FALSE.," << endl;
		ofs << "  wfdt   = 2.0D0," << endl;
		ofs << "  wf_q   = 500.D0," << endl;
		ofs << "  wf_friction = 0.3d0," << endl;

		ofs << "ATOMIC_SPECIES" << endl;
		ofs << "O   15.9994    O_HSCV_PBE-1.0.UPF" << endl;
		ofs << "H   2.01410178 H_HSCV_PBE-1.0.UPF" << endl;

		ofs << "ATOMIC_POSITIONS {bohr}" << endl;
		ifs >> snapshot_index >> snapshot_time;
		ofs << "! " << snapshot_index << " " << snapshot_time << endl;
		ofs << setprecision(14);
		for(int ia = 0; ia <INPUT.natom1; ia++)
		{
			ifs >> x >> y >> z;
			ofs << "O    " << x << "    " << y << "    " << z << endl;
		}
		for(int ia = 0; ia <INPUT.natom2; ia++)
		{
			ifs >> x >> y >> z;
			ofs << "H    " << x << "    " << y << "    " << z << endl;
		}
		ofs.close();
		cout << "igeo = " << igeo << endl;
	}
	ifs.close();
}
