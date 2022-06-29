#include "fp_check.h"
#include "PT_snapshot.h"

fp_check::fp_check()
{}

fp_check::~fp_check()
{}

void fp_check::Routine()
{
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
			
		}
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) 
		{
			cel.clean();//qianrui add in 2020-1-7
			continue;
		}
		cout << "snapshot " << igeo << endl;
		string iigeo = to_string(igeo);
		ofstream ofs_scf(iigeo + ".in");
		ofstream ofs_json(iigeo + ".json");
		//print_input(cel, ofs_scf, ofs_json, iigeo);
		print_input_cp(cel, ofs_scf, ofs_json, iigeo);

	}
}



void fp_check::print_input(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index)
{
	ofs_scf << setprecision(14);
	ofs_scf << "&CONTROL !" << cel.snapshot_index << " " << cel.snapshot_time << endl;
	ofs_scf << "title = \"water 63mol+OH- SCAN 330K\"," << endl;
	ofs_scf << "calculation = 'scf'," << endl;
	ofs_scf << "restart_mode = \"from_scratch\"," << endl;
	ofs_scf << "prefix = \"oh-_"+out_index+"\"," << endl;
	ofs_scf << "pseudo_dir = '"+INPUT.pseudo_in+"'," << endl;
	ofs_scf << "outdir = './SimulationFiles/'," << endl;
	ofs_scf << "nstep = 1," << endl;
	ofs_scf << "isave = 1," << endl;
	ofs_scf << "dt = 20.D0," << endl;
	ofs_scf << "tstress = .true.," << endl;
	ofs_scf << "tprnfor = .true.," << endl;
	ofs_scf << "iprint = 1," << endl;
	ofs_scf << "etot_conv_thr = 1.D-6," << endl;
	ofs_scf << "max_seconds = 43000," << endl;
	ofs_scf << "/" << endl;
	ofs_scf << "&SYSTEM" << endl;
	ofs_scf << "ibrav     = 1," << endl;
	ofs_scf << "celldm(1) = 23.5170," << endl;
	ofs_scf << "nat       = 191," << endl;
	ofs_scf << "ntyp      = 2," << endl;
	ofs_scf << "ecutwfc   = 85.00D0," << endl;
	ofs_scf << "ecfixed   = 130.D0, ! effective cutoff value" << endl;
	ofs_scf << "q2sigma   = 15.D0," << endl;
	ofs_scf << "qcutz     = 200.D0," << endl;
	ofs_scf << "input_dft = 'scan'," << endl;
	ofs_scf << "tot_charge    = -1," << endl;
	ofs_scf << "/" << endl;
	ofs_scf << "&ELECTRONS" << endl;
	ofs_scf << "ortho_max             = 800," << endl;
	ofs_scf << "ortho_eps             = 1.D-7," << endl;
	ofs_scf << "electron_maxstep      = 500," << endl;
	ofs_scf << "/" << endl;
	ofs_scf << "&IONS" << endl;
	ofs_scf << "tempw            = 330," << endl;
	ofs_scf << "/" << endl;
	ofs_scf << "ATOMIC_SPECIES" << endl;
	ofs_scf << "O   15.9994    O_HSCV_PBE-1.0.UPF" << endl;
	ofs_scf << "H   2.01410178 H_HSCV_PBE-1.0.UPF" << endl;
	ofs_scf << "ATOMIC_POSITIONS {bohr}" << endl;

	ofs_geo_json << setprecision(14);
	ofs_geo_json << "[";

	for(int it = 0; it < INPUT.ntype; it++)
	{
		for(int ia = 0; ia < cel.atom[it].na; ia++)
		{
			ofs_scf << cel.atom[it].id << " "  << cel.atom[it].pos[ia].x/BOHR << " " << cel.atom[it].pos[ia].y/BOHR << " " << cel.atom[it].pos[ia].z/BOHR << endl;

			ofs_geo_json << "[" << cel.atom[it].pos[ia].x << ", " << cel.atom[it].pos[ia].y << ", " << cel.atom[it].pos[ia].z << "]";
			if(ia != cel.atom[it].na-1 or it != INPUT.ntype-1) ofs_geo_json << ", " << endl;
		}
	}
	ofs_geo_json << "]" << endl;
	ofs_scf.close();
	ofs_geo_json.close();
}

void fp_check::print_input_cp(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index)
{
	ofs_scf << setprecision(14);
	ofs_scf << "&CONTROL !" << cel.snapshot_index << " " << cel.snapshot_time << endl;
	ofs_scf << "title = \"water 63mol+OH- SCAN 330K\"," << endl;
	ofs_scf << "calculation = 'cp'," << endl;
	ofs_scf << "restart_mode = \"from_scratch\"," << endl;
	ofs_scf << "prefix = \"oh-_"+out_index+"\"," << endl;
	ofs_scf << "pseudo_dir = '"+INPUT.pseudo_in+"'," << endl;
	ofs_scf << "outdir = './SimulationFiles/'," << endl;
	ofs_scf << "nstep = 100000," << endl;
	ofs_scf << "isave = 100," << endl;
	ofs_scf << "dt = 2.D0," << endl;
	ofs_scf << "tstress = .true.," << endl;
	ofs_scf << "tprnfor = .true.," << endl;
	ofs_scf << "iprint = 1," << endl;
	ofs_scf << "etot_conv_thr = 1.D-6," << endl;
	ofs_scf << "/" << endl;
	
	ofs_scf << "&SYSTEM" << endl;
	ofs_scf << "ibrav     = 1," << endl;
	ofs_scf << "celldm(1) = 23.5170," << endl;
	ofs_scf << "nat       = 191," << endl;
	ofs_scf << "ntyp      = 2," << endl;
	ofs_scf << "ecutwfc   = 85.00D0," << endl;
	ofs_scf << "ecfixed   = 130.D0, ! effective cutoff value" << endl;
	ofs_scf << "q2sigma   = 15.D0," << endl;
	ofs_scf << "qcutz     = 200.D0," << endl;
	ofs_scf << "input_dft = 'scan'," << endl;
	ofs_scf << "tot_charge    = -1," << endl;

	ofs_scf << "/" << endl;
	ofs_scf << "&ELECTRONS" << endl;
	ofs_scf << "emass                 = 100.D0" << endl;
	ofs_scf << "emass_cutoff          = 25.D0" << endl;
	ofs_scf << "electron_dynamics     = 'damp'" << endl;
	ofs_scf << "electron_damping      = 0.20" << endl;
	ofs_scf << "ortho_max             = 800," << endl;
	ofs_scf << "ortho_eps             = 1.D-7," << endl;

	ofs_scf << "/" << endl;
	ofs_scf << "&IONS" << endl;
	ofs_scf << "ion_dynamics     = \"none\"," << endl;
	ofs_scf << "ion_temperature  = 'not_controlled'," << endl;
	ofs_scf << "ion_radius(1) = 1.4" << endl;
	ofs_scf << "ion_radius(2) = 1.4" << endl;

	ofs_scf << "/" << endl;
	ofs_scf << "ATOMIC_SPECIES" << endl;
	ofs_scf << "O   15.9994    O_HSCV_PBE-1.0.UPF" << endl;
	ofs_scf << "H   2.01410178 H_HSCV_PBE-1.0.UPF" << endl;
	ofs_scf << "ATOMIC_POSITIONS {bohr}" << endl;

	ofs_geo_json << setprecision(14);
	ofs_geo_json << "[";

	for(int it = 0; it < INPUT.ntype; it++)
	{
		for(int ia = 0; ia < cel.atom[it].na; ia++)
		{
			ofs_scf << cel.atom[it].id << " "  << cel.atom[it].pos[ia].x/BOHR << " " << cel.atom[it].pos[ia].y/BOHR << " " << cel.atom[it].pos[ia].z/BOHR << endl;

			ofs_geo_json << "[" << cel.atom[it].pos[ia].x << ", " << cel.atom[it].pos[ia].y << ", " << cel.atom[it].pos[ia].z << "]";
			if(ia != cel.atom[it].na-1 or it != INPUT.ntype-1) ofs_geo_json << ", " << endl;
		}
	}
	ofs_geo_json << "]" << endl;
	ofs_scf.close();
	ofs_geo_json.close();
}
