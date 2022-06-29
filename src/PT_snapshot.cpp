// The purpose of this subroutine is to find out snapshots going through PTs and
// transform them into QE-scf input scripts & json list files. renxi 20200613
// nPT in INPUT denotes # of PT needed to be converted
#include "PT_snapshot.h"
PT_snapshot::PT_snapshot()
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
	snapshot_printed = new int[INPUT.nPT];
	snapshot_printed_pr = new int[INPUT.nPT];
	snapshot_printed_aft = new int[INPUT.nPT];
}

PT_snapshot::~PT_snapshot()
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
	delete[] snapshot_printed;
	delete[] snapshot_printed_pr;
	delete[] snapshot_printed_aft;
}

void PT_snapshot::Routine()
{
	ifstream ifs_PT("mj.dat");
	string useless;
	double BOHR = 0.52917721092;
	ofstream ofs_H("H_index.dat");
	ofstream ofs_O("O_index.dat");
	ofstream ofs_time("time.dat");
	cout << "INPUT.geo_target_aft " << INPUT.geo_target_aft << endl;
	for(int i = 0; i < INPUT.nPT+1; i++)
	{
		ifs_PT >> index[i] >> PT_start_snapshot_index[i] >> PT_start_snapshot_time[i] >> ion_index[i] >> nacc_start[i] >> ndon_start[i] >> PT_end_snapshot_index[i] >> PT_end_snapshot_time[i] >> useless >> nacc_end[i] >> ndon_end[i] >> H_index[i] >> useless;
		cout << i << " " << PT_start_snapshot_index[i] << " " << PT_end_snapshot_index[i] << endl;
		assert(index[i]==i+1);
		snapshot_printed[i] = 0;
		snapshot_printed_pr[i] = 0;
		snapshot_printed_aft[i] = 0;
		ofs_H << H_index[i] << endl;
		ofs_O << ion_index[i] << endl;
	}
	ofs_H.close();
	ofs_O.close();
	ifs_PT.close();
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
			cout << "snapshot " << igeo << endl;
			continue;
		}
		cout << "snapshot " << igeo << endl;
		cout << "snapshot_index " << cel.snapshot_index << endl;
		for(int iPT = 0; iPT < INPUT.nPT; iPT++)
		{
			if(cel.snapshot_index >= PT_end_snapshot_index[iPT] - INPUT.geo_target and cel.snapshot_index < PT_end_snapshot_index[iPT] and cel.snapshot_index%INPUT.nbin == 0)
			{
				snapshot_printed_pr[iPT]++;
				string out_index_pr = "-"+to_string(iPT+1)+"_"+to_string(snapshot_printed_pr[iPT]);
				ofs_time << out_index_pr << " " << cel.snapshot_index << " " << cel.snapshot_time << endl;
				ofstream ofs_scf_pre("oh-_" + out_index_pr + ".in");
				ofstream ofs_geo_json_pre("oh-_geo_" + out_index_pr + ".json");
				print_input(cel, ofs_scf_pre, ofs_geo_json_pre, out_index_pr);
			}
			if(cel.snapshot_index >= PT_end_snapshot_index[iPT] and cel.snapshot_index <= PT_start_snapshot_index[iPT+1])
			{
				snapshot_printed[iPT]++;
				string out_index = to_string(iPT+1) + "_" + to_string(snapshot_printed[iPT]);
				ofs_time << out_index << " " << cel.snapshot_index << " " << cel.snapshot_time << endl;
				ofstream ofs_scf("oh-_" + out_index + ".in");
				ofstream ofs_geo_json("oh-_geo_" + out_index + ".json");
				print_input(cel, ofs_scf, ofs_geo_json, out_index);
			}
			if(cel.snapshot_index <= PT_start_snapshot_index[iPT+1] + INPUT.geo_target_aft and cel.snapshot_index > PT_start_snapshot_index[iPT+1] and cel.snapshot_index%INPUT.nbin == 0)
			{
				snapshot_printed_aft[iPT]++;
				string out_index_aft = "+"+to_string(iPT+1)+"_"+to_string(snapshot_printed_aft[iPT]);
				ofs_time << out_index_aft << " " << cel.snapshot_index << " " << cel.snapshot_time << endl;
				ofstream ofs_scf_aft("oh-_" + out_index_aft + ".in");
				ofstream ofs_geo_json_aft("oh-_geo_" + out_index_aft + ".json");
				print_input(cel, ofs_scf_aft, ofs_geo_json_aft, out_index_aft);
			}

		}
	}
	ofs_time.close();
}

void PT_snapshot::print_input(const Cell &cel, ofstream &ofs_scf, ofstream &ofs_geo_json, string &out_index)
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
