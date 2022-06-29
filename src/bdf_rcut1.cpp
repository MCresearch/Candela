#include "bdf_rcut1.h"
#include "gfun.h"

BDF_rcut1::BDF_rcut1(){}

BDF_rcut1::~BDF_rcut1(){}

void BDF_rcut1::Routine()
{
    TITLE("BDF_rcut1","Routine");

    this->bdf_angle = new double[int(180/INPUT.bdf_dtheta)+1];
    for (int iangle=0; iangle < int(180/INPUT.bdf_dtheta)+1; iangle++)
    {
        bdf_angle[iangle] = 0;
    }

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
        this->calc(cel);
        cel.clean();
	}//igeo
    this->output();
    return;
}


void BDF_rcut1::calc(CellFile &cel)
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
    
    Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;
    HBs::setup_water(cel, water);

    for (int iwater0=0; iwater0<cel.atom[ito].na; iwater0++)
    {
        for (int iwater1=0; iwater1<cel.atom[ito].na; iwater1++)
        {
            if (iwater0 == iwater1) continue;
            if (distance(cel.atom[ito].pos[iwater0], cel.atom[ito].pos[iwater1], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3) > INPUT.bdf_rcut) continue;

            for (int iwater2=0; iwater2<cel.atom[ito].na; iwater2++)
            {
                if (iwater0 == iwater2 or iwater1 == iwater2) continue;
                if (distance(cel.atom[ito].pos[iwater0], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3) > INPUT.bdf_rcut) continue;
            
                double bond_angle = HBs::angle(cel, cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater0], cel.atom[ito].pos[iwater2]);

                int which = int(bond_angle/INPUT.bdf_dtheta);

                if (which < int(180/INPUT.bdf_dtheta)+1)
                {
                    this->bdf_angle[which]++;
                }
            }
        }
    }
    delete[] water;
}

void BDF_rcut1::output()
{
    double summ = 0;
    for (int iangle=0; iangle<int(180/INPUT.bdf_dtheta)+1; iangle++)
    {
        summ += this->bdf_angle[iangle];
    }
    ofstream ofs("bdf_angle.txt");
    for (int iangle=0; iangle<int(180/INPUT.bdf_dtheta)+1; iangle++)
    {
        ofs << iangle*INPUT.bdf_dtheta << " " << this->bdf_angle[iangle]/summ/INPUT.bdf_dtheta << endl;
    }
    ofs.close();
    return;
}