#include "cellFile.h"
#include "input.h"
#include "math.h"
#include "tune_stru.h"
#include "HBs.h"

Tune_Stru::Tune_Stru()
{
}

Tune_Stru::~Tune_Stru()
{
}

void Tune_Stru::Routine()
{
	TITLE("Tune_Stru","Routine");
	
	cout << "Tune_Stru structure" << endl;

	ofstream ofs("constrain.dat");
	ofstream ofs2("tuned_stru.xyz");

	this->count_geometry_number=0;
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

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "igeo=" << igeo << endl;

		adjust_structure(cel, igeo, ofs, ofs2);
	}	
	
	ofs2.close();

	return;
}


void Tune_Stru::adjust_structure(const Cell &cel, const int &igeo, ofstream &ofs, ofstream &ofs2)
{
    int ito=-1;
    int ith=-1;
    int itc=-1;
    for(int it=0;it <INPUT.ntype; ++it)
    {
        if(cel.atom[it].id=="O") ito=it;
        else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
        else if(cel.atom[it].id=="C") itc=it;
    }
    if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
    if(INPUT.ntype==3){ assert(itc>=0); }

    Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;

    HBs::setup_water(cel, water);

	bool* labelH = new bool[cel.atom[ith].na];
	for(int i=0; i<cel.atom[ith].na; ++i)
	{
		labelH[i]=false;
	}

	ofs << "CONSTRAINTS" << endl;
	ofs << "number_of_constrains" << endl;	

	ofs2 << INPUT.natom << endl;
	ofs2 << "title" << endl;

	int iat=0;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		++iat;
		int io=iat;

		double x = cel.atom[ito].pos[ia].x;
		double y = cel.atom[ito].pos[ia].y;
		double z = cel.atom[ito].pos[ia].z;

		assert(INPUT.factor>0);
		double newX = x; 
		double newY = y;
		double newZ = z * INPUT.factor;

		ofs2 << cel.atom[ito].id << " " << newX << " " << newY << " " << newZ << endl;

		for(int ia2=0; ia2<water[ia].nH; ++ia2)
		{
			++iat;

			const int indexH = water[ia].indexH[ia2];  
			labelH[indexH]=true;
	
			double tmpX = cel.atom[ith].pos[indexH].x; 
			double tmpY = cel.atom[ith].pos[indexH].y; 
			double tmpZ = cel.atom[ith].pos[indexH].z; 
	
			double dx = shortest(tmpX,x,INPUT.celldm1);
			double dy = shortest(tmpY,y,INPUT.celldm2);
			double dz = shortest(tmpZ,z,INPUT.celldm3);

			double Hx = dx + newX;			
			double Hy = dy + newY;			
			double Hz = dz + newZ;	

			ofs2 << cel.atom[ith].id << " " << Hx << " " << Hy << " " << Hz << endl; // .xyz		

			// difficult to understand but less coding efforts
			if(water[ia].nH==2 and ia2==1) ofs << " " << iat << endl;
			ofs << "'distance' " << io << " " << iat << endl;
			if(water[ia].nH==2 and ia2==0) ofs << "'planar_angle' " << iat << " " << io; 
			//if(water[ia].nH==3) ofs << "#########" << endl;

		}
	}

	// print out carbon information
	for(int ia3=0; ia3<cel.atom[itc].na; ++ia3)
	{
		double Cx = cel.atom[itc].pos[ia3].x; 
		double Cy = cel.atom[itc].pos[ia3].y; 
		double Cz = cel.atom[itc].pos[ia3].z; 

		ofs2 << cel.atom[itc].id 
			<< " " << Cx  
			<< " " << Cy 
			<< " " << Cz * INPUT.factor << endl; 
	}

	// print out rest hydrogen information
	for(int ia4=0; ia4<cel.atom[ith].na; ++ia4)
	{
		if( labelH[ia4]==false )
		{
			ofs2 << cel.atom[ith].id 
				<< " " << cel.atom[ith].pos[ia4].x 
				<< " " << cel.atom[ith].pos[ia4].y 
				<< " " << cel.atom[ith].pos[ia4].z * INPUT.factor 
				<< endl;
		}
	}



	delete[] water;
	delete[] labelH;
	
	return;
}
