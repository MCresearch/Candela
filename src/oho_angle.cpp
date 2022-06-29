#include "oho_angle.h"
oho_angle::oho_angle(){}

oho_angle::~oho_angle(){}

void oho_angle::Routine()
{

    assert(INPUT.dtheta > 0);
    assert(INPUT.system == "water");
    assert(INPUT.rcut > 0);
    this->ntheta = int(180/INPUT.dtheta);
    this->theta_distr = new double[ntheta];
    for (int ia=0; ia<ntheta; ia++) theta_distr[ia] = 0;

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
			cel.clean();//qianrui add in 2020-1-7
			continue;
		}
		++count_geometry_number;
		cout << "snapshot " << igeo << endl;

		int ito = -1;
		int ith = -1;
		for(int it=0; it<INPUT.ntype; it++)
		{
			if(cel.atom[it].id=="O"){ito = it;}
			else if(cel.atom[it].id=="H" or cel.atom[it].id=="D"){ith = it;}
		}
		assert(ito>=0 and ith>=0);
		Water *water = new Water[cel.atom[ito].na];
		Water::nions = 0;
		HBs::setup_water(cel, water);

        calc(cel, water, ito, ith);
    }
    normalize();
    ofstream ofs("oho_angle.txt");
    for (int ia=0; ia<this->ntheta; ia++)
    {
        ofs << setprecision(12) << ia*INPUT.dtheta << " " << this->theta_distr[ia] << endl;
    }
    ofs.close();
}

void oho_angle::calc(CellFile &cel, Water* water, int &ito, int &ith)
{
    for (int ia=0; ia<cel.atom[ito].na; ia++)
    {
        for (int ia2=0; ia2<cel.atom[ito].na; ia2++)
        {
            if (ia == ia2) continue;
            if (distance(cel.atom[ito].pos[ia], cel.atom[ito].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3) > INPUT.rcut) continue;
            double angle1 = HBs::angle(cel, cel.atom[ito].pos[ia], cel.atom[ith].pos[water[ia].indexH[0]], cel.atom[ito].pos[ia2]);
            //cout << ia << " " << ia2 << " " << 0 << " " << angle1 << endl;
            this->theta_distr[int(angle1/INPUT.dtheta)]++;
            angle1 = HBs::angle(cel, cel.atom[ito].pos[ia], cel.atom[ith].pos[water[ia].indexH[1]], cel.atom[ito].pos[ia2]);
            //cout << ia << " " << ia2 << " " << 1 << " " << angle1 << endl;
            this->theta_distr[int(angle1/INPUT.dtheta)]++;
        }
    }
    return;
}

void oho_angle::normalize()
{
    double sum = 0;
    for (int ia=0; ia<this->ntheta; ia++)
    {
        sum += this->theta_distr[ia];
    }
    sum *= INPUT.dtheta;
    for (int ia=0; ia<this->ntheta; ia++)
    {
        this->theta_distr[ia] /= sum;
    }
    return;
}