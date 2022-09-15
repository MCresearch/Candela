#include "cellFile.h"
#include "input.h"
#include "write.h"
#include "stdlib.h"

void Write::Routine()
{

	cout<<"Writing GEO input file..."<<endl;
	assert(INPUT.geo_interval>0);

	cout << INPUT.geo_1 << " " << INPUT.geo_2 << endl;
	int count=1;

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo) 
	{
		CellFile cel;

		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		if( !CellFile::ReadGeometry( cel ) ) continue;
		if(cel.read_and_used==false) 
		{
			cel.clean();
			continue;
		}
		ifstream ifs(INPUT.headfile.c_str());
		string txt;
		string countstr=int2str(count);
		string command="test -e "+countstr+" ||mkdir "+countstr;
		cout<<command<<endl;
		system(command.c_str());
		string filename="./"+countstr+"/"+INPUT.geo_out;
		ofstream ofs(filename.c_str());
		while(getline(ifs,txt))
		{
			if(ifs.eof()) break;
			ofs<<txt<<endl;
		}
		ifs.close();
		//cout<<INPUT.ntype<<' '<<cel.atom[0].na<<endl;
		if(!INPUT.write_cartesian)
			cout<<cel.a1.x<<' '<<cel.a2.y<<' '<<cel.a3.z<<endl;//test
		for(int it=0;it<INPUT.ntype;it++)
		{
			for(int ia=0;ia<cel.atom[it].na;ia++)
			{
				if(INPUT.write_cartesian)
				ofs<<cel.atom[it].id<<'\t'<<cel.atom[it].pos[ia].x*INPUT.factor<<'\t'<<cel.atom[it].pos[ia].y*INPUT.factor<<'\t'<<cel.atom[it].pos[ia].z*INPUT.factor<<"\t1\t1\t1"<<endl;
				else
				ofs<<cel.atom[it].id<<'\t'<<cel.atom[it].pos[ia].x/cel.a1.x<<'\t'<<cel.atom[it].pos[ia].y/cel.a2.y<<'\t'<<cel.atom[it].pos[ia].z/cel.a3.z<<"\t1\t1\t1"<<endl;
				
			}
		}
		ifs.open(INPUT.tailfile.c_str());
		while(getline(ifs,txt))
		{
			if(ifs.eof()) break;
			ofs<<txt<<endl;
		}
		ifs.close();
		ofs.close();
        cel.clean();
		count++;
	}
}
