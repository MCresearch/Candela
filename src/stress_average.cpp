#include "stress_average.h"
#include "input.h"
#include "gfun.h"
stress_average::stress_average(){}

stress_average::~stress_average(){}

void stress_average::Routine()
{
	int ngeo = INPUT.geo_2-INPUT.geo_1+1;
	stress_list = new double[ngeo+10]();
	stringstream ss;
	ifstream ifs;
	ss << INPUT.geo_directory;
	cout << " Read stress : " << ss.str() << endl;
	ifs.open(ss.str().c_str());
	if(!ifs)
	{
		cout << "could not find the .str file: " << INPUT.geo_directory << endl;
		exit(0);
	} 
	cout << " File name is " << ss.str() << endl;
	//cout << " CellFile::file_open = " << file_open << endl;
	ofstream ofs("stress.txt");
	int istress = 0;
	double snapshot_time;
	int snapshot_index;
	string useless;
	double str_x, str_y, str_z;
	double stress;
	double average_stress = 0;
	double standard_error = 0;
	for(int igeo = INPUT.geo_1; igeo<=INPUT.geo_2; igeo++)
	{
		ifs >> snapshot_index >> snapshot_time;
		ifs >> str_x >> useless >> useless;
		ifs >> useless >> str_y >> useless;
		ifs >> useless >> useless >> str_z;
		stress = (str_x+str_y+str_z)/3;
		stress_list[istress] = stress;
		cout << "igeo = " << igeo << endl;
		cout << "snapshot_index " << snapshot_index << endl;
		cout << "snapshot_time " << snapshot_time << endl; 
		ofs << igeo << " " << snapshot_time << " " << stress << endl;
		average_stress+=stress;
		istress++;
	}
	average_stress /= istress;
	for(int is = 0; is<istress; is++)
	{
		standard_error += pow(stress_list[is], 2);
	}
	standard_error = sqrt(standard_error/istress - pow(average_stress, 2));
	cout << "average_stress " << average_stress << endl;
	cout  << "standard_error " << standard_error << endl;
	ofs.close();
	ifs.close();
}
