#include "input.h"
#include "cellFile.h"
#include "band_gap.h"

band_gap::band_gap()
{}

band_gap::~band_gap()
{}

void band_gap::Routine()
{
	TITLE("band_gap", "Routine");

	ofs_running << "calculate and output the band gap of every snapshot." << endl;
	int HO = INPUT.HO;
	int LU = INPUT.LU;
	ofstream ofs("band_gap.txt");
	int count_geometry_number = 0;
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
			cel.clean(); // renxi added 20200614
			continue;
		}
		++count_geometry_number;
		cout << "snapshot " << igeo << endl;
		//cout << cel.eig[0] << " " << cel.eig[289] << endl;
		double band_gap = cel.eig[LU]-cel.eig[HO];
		ofs << cel.snapshot_index << " " << cel.snapshot_time << " " << cel.eig[HO] << " " << cel.eig[LU] << " " << band_gap << endl;
	}
	ofs.close();
}
