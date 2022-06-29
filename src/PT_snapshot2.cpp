#include "PT_snapshot2.h"
#include "PT_snapshot.h"
// Usage of this subroutine:
// output snapshots around PT in pw.x's format and json format. renxi 20200622
// INPUT.output_interval: length of time interval around PT where snapshots are printed;
// INPUT.nsnapshot: number of snapshots printed every PT;
// INPUT.nPT: number of PT we want to deal with.

PT_snapshot2::PT_snapshot2()
{
	index = new int[INPUT.nPT];
	PT_start_snapshot_index = new int[INPUT.nPT];
	PT_start_snapshot_time = new double[INPUT.nPT];
	PT_end_snapshot_time = new double[INPUT.nPT];
	PT_end_snapshot_index = new int[INPUT.nPT];
	ion_index = new int[INPUT.nPT];
	nacc_start = new int[INPUT.nPT];
	ndon_start = new int[INPUT.nPT];
	nacc_end = new int[INPUT.nPT];
	ndon_end = new int[INPUT.nPT];
	H_index = new int[INPUT.nPT];
}

PT_snapshot2::~PT_snapshot2()
{
	delete[] index;
	delete[] PT_start_snapshot_index; // the start here acturally does not mean start of PT, rather, start of rest
	delete[] PT_start_snapshot_time;
	delete[] PT_end_snapshot_time;
	delete[] PT_end_snapshot_index;
	delete[] ion_index;
	delete[] nacc_start;
	delete[] ndon_start;
	delete[] nacc_end;
	delete[] ndon_end;
	delete[] H_index;
}

void PT_snapshot2::Routine()
{
	ifstream ifs_PT("mj.dat");
	string useless;
	double BOHR = 0.52917721092;
	double snapshot_time;
	int snapshot_index;
	ofstream ofs("print_record.dat");
	ofs << "index snapshot_index PT_ss_index" << endl;
	int iindex = 0;
	double tot_energy;
	for(int i = 0; i < INPUT.nPT+1; i++)
	{
		ifs_PT >> index[i] >> PT_start_snapshot_index[i] >> PT_start_snapshot_time[i] >> ion_index[i] >> nacc_start[i] >> ndon_start[i] >> PT_end_snapshot_index[i] >> PT_end_snapshot_time[i] >> useless >> nacc_end[i] >> ndon_end[i] >> H_index[i] >> useless;
		cout << i << " " << PT_start_snapshot_index[i] << " " << PT_end_snapshot_index[i] << endl;
		assert(index[i]==i+1);
	}
	ifs_PT.close();
	// INPUT.output_interval: length of time interval around PT where snapshots are printed;
	// INPUT.nsnapshot: number of snapshots printed every PT;
	// INPUT.nPT: number of PT we want to deal with.

	double ss_time_interval = INPUT.output_interval/INPUT.nsnapshot;
	int ss_index_interval = ss_time_interval/INPUT.msd_dt;

	if (INPUT.func ==1 )
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
				for (int iPT=0; iPT <= INPUT.nPT; iPT++)
				{
					if((cel.snapshot_index-PT_end_snapshot_index[iPT])%ss_index_interval==0 and 
						(abs((cel.snapshot_index-PT_end_snapshot_index[iPT])*INPUT.msd_dt) <= INPUT.output_interval/2))
					{
						ofs << iindex << " " << cel.snapshot_index << " " << PT_end_snapshot_index[iPT] << endl;
						string iindexs = to_string(iindex);
						ofstream ofs_scf(iindexs + ".in");
						ofstream ofs_geo(iindexs + ".json");
						//print_input(cel, ofs_scf, ofs_geo, iindexs);
						print_input(cel, ofs_scf, ofs_geo, iindexs);
						iindex++;
					}
				}
		}
	}	
	if (INPUT.func == 2)
	{
		ifstream ifs_eng(INPUT.geo_directory);
		ofstream ofs_eng("PT_tot_energy.dat");
		for (int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
		{

			ifs_eng >> snapshot_index >> snapshot_time >> useless >> useless >> useless >> tot_energy >> useless >> useless >> useless >> useless >> useless;
			for (int iPT=0; iPT <= INPUT.nPT; iPT++)
			{
				if((snapshot_index-PT_end_snapshot_index[iPT])%ss_index_interval==0 and 
					(abs(snapshot_time-PT_end_snapshot_time[iPT]) <= INPUT.output_interval/2))
				{
					ofs_eng << setprecision(14) << iPT << " " << snapshot_time-PT_end_snapshot_time[iPT] << " " << tot_energy << endl;
				}
			}
			cout << "geo = " << igeo << endl;
		}
		ifs_eng.close();
		ofs_eng.close();
	}
	if (INPUT.func == 3)
	{
		ifstream ifs_eng(INPUT.geo_directory);
		ofstream ofs_eng("PT_tot_energy.dat");
		while (ifs_eng.good())
		{
			ifs_eng >> useless;
			if (useless == "Volume")
			{
				break;
			}
		}
		for (int igeo = INPUT.geo_1; igeo <= INPUT.geo_2; ++igeo)
		{
			ifs_eng >> snapshot_index >> tot_energy >> useless >> useless >> useless >> useless >> useless;
			for (int iPT=0; iPT <= INPUT.nPT; iPT++)
			{
				if (abs((snapshot_index-PT_end_snapshot_index[iPT])*INPUT.msd_dt) <= INPUT.output_interval/2)
				{
					ofs_eng << setprecision(20) << tot_energy << endl;
				}
			}
			cout << "geo = " << igeo << endl;
		}
		ifs_eng.close();
		ofs_eng.close();
	}
	ofs.close();
}

void PT_snapshot2::print_input(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index)
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

void PT_snapshot2::print_input_cp(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index)
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
