#include "reorganize.h"
#include "input.h"
#include "cellFile.h"
#include "HBs.h"

Reorganize::Reorganize() 
{
}

Reorganize::~Reorganize() 
{
}

void Reorganize::Routine()
{
	cout << " Reorganize the snapshots" << endl;

	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

	ofstream ofsp("new.pos");
	ofstream ofsw("new.wfc");
	ofstream ofsc("new.cel");
	ofstream ofsx("new.xyz");
	ofsp << setprecision(12);
	ofsw << setprecision(12);
	ofsc << setprecision(12);
	ofsx << setprecision(6);

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//ofs_running << "geometry " << igeo%INPUT.geo_interval << endl;
		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;


		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) )
		{
			cel.clean();
			continue;
		}
		if(cel.read_and_used==false) 
		{
			cel.clean();
			continue;
		}
		++count_geometry_number;
		cout << "snapshot " << igeo << endl;

		print_out(cel, ofsp, ofsw, ofsc, ofsx);
		cel.clean();
	}

	ofsp.close();
	ofsw.close();
	ofsc.close();
	ofsx.close();
	return;
}

void Reorganize::print_out(Cell &cel, ofstream &ofsp, ofstream &ofsw, ofstream &ofsc, ofstream &ofsx)
{
	if(cel.snapshot_index%INPUT.nbin!=0) return;

	const double bohr=0.529177;

	ofsp << cel.snapshot_index << " " << cel.snapshot_time << endl; 
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			ofsp << cel.atom[it].pos[ia].x/bohr 
				<< " " << cel.atom[it].pos[ia].y/bohr 
				<< " " << cel.atom[it].pos[ia].z/bohr << endl;
		}
	}

	// ofsx: xyz file
	int nat=0;
	for(int it=0; it<INPUT.ntype; ++it) nat+=cel.atom[it].na;
	ofsx << nat << endl;
	ofsx << cel.snapshot_index << " " << cel.snapshot_time << endl; 
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			ofsx << cel.atom[it].id
				<< " " << cel.atom[it].pos[ia].x 
				<< " " << cel.atom[it].pos[ia].y 
				<< " " << cel.atom[it].pos[ia].z << endl;
		}
	}

	ofsw << cel.snapshot_index << " " << cel.snapshot_time << endl; 
	for(int ib=0; ib<INPUT.nbands; ++ib)
	{
		// distance between Wannier centers and oxygen atoms
		ofsw << cel.wan_centers[ib].x/bohr 
			<< " " << cel.wan_centers[ib].y/bohr 
			<< " " << cel.wan_centers[ib].z/bohr << endl;
	}

	ofsc << cel.snapshot_index << " " << cel.snapshot_time << endl;
	ofsc << cel.a1.norm()/bohr << " 0 0" << endl;
	ofsc << "0 " << cel.a2.norm()/bohr << " 0" << endl;
	ofsc << "0 0 " << cel.a3.norm()/bohr << endl;

	return;
} 
