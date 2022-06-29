#include "Wan_centers_stat.h"
#include "input.h"

Wan_centers_stat::Wan_centers_stat(){}

Wan_centers_stat::~Wan_centers_stat(){}

void Wan_centers_stat::Routine()
{

    if (INPUT.wannier_file == "none")
    {
        cout << "Wannier file is a must for this analysis." << endl;
        exit(0);
    }

	this->nangle = (int) (INPUT.theta/INPUT.dtheta)+1;
	this->nr = (int) (INPUT.rcut/INPUT.dr)+1;
	this->lone_angle_distr = new double[nangle];
	for(int ia = 0; ia < nangle; ia++)
	{
		lone_angle_distr[ia] = 0.0;
	}
	this->lone_2D_distr = new double*[nr];
	for (int ir=0; ir<nr; ir++)
	{
		lone_2D_distr[ir] = new double[nr];
		for (int ir2=0; ir2<nr; ir2++)
		{
			lone_2D_distr[ir][ir2] = 0;
		}
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

        this->Wan_angle_distr(cel);
	}// igeo

	out(count_geometry_number);

	return;
}

void Wan_centers_stat::Wan_angle_distr(Cell &cel)
// This function calculates the angle distribution of the two lone pairs and their 2D distribution
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
    int** bond_wan_index = new int*[cel.atom[ito].na];
    int** lone_wan_index = new int*[cel.atom[ito].na];
    for (int ia=0; ia < cel.atom[ito].na; ia++)
    {
        bond_wan_index[ia] = new int[3];
        lone_wan_index[ia] = new int[3];
    }

	this->allocate_wan(water, cel, bond_wan_index, lone_wan_index, ito, ith);
	this->calculate(water, cel, bond_wan_index, lone_wan_index, ito, ith);

	delete[] water;
    for (int ia=0; ia < cel.atom[ito].na; ia++)
    {
        delete[] bond_wan_index[ia];
        delete[] lone_wan_index[ia];
    }

    delete[] bond_wan_index;
    delete[] lone_wan_index;
    return;
}

void Wan_centers_stat::allocate_wan(const Water* water, const Cell &cel, int** &bond_wan_index, int** &lone_wan_index, const int &ito, const int &ith)
{
	for(int iwater = 0; iwater<cel.atom[ito].na; iwater++)
	{
		int nwan = 0; 
		int wan_index[4];
		for(int iwan = 0; iwan < INPUT.nbands; iwan++)
		{
			if (distance(cel.atom[ito].pos[iwater], cel.wan_centers[iwan], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3) < 1)
			{
				wan_index[nwan] = iwan;
				nwan++;
				if (nwan > 4)
				{
					cout << "nwan cannot be bigger than 4. Check." << endl;
					exit(0);
				}
			}
		}
		assert(nwan == 4);
		int nlone = 0;
		int nbond = 0;
		for (int iwan = 0; iwan < 4; iwan++)
		{
			double dis0 = distance(cel.atom[ith].pos[water[iwater].indexH[0]], cel.wan_centers[wan_index[iwan]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
			double dis1 = distance(cel.atom[ith].pos[water[iwater].indexH[1]], cel.wan_centers[wan_index[iwan]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
			if (dis0 < 0.9 or dis1 < 0.9)
			{
				bond_wan_index[iwater][nbond] = wan_index[iwan];
				nbond++;
			}
			else
			{
				lone_wan_index[iwater][nlone] = wan_index[iwan];
				nlone++;
			}
		}//iwan
		if (nlone != 2)
		{
			cout << "nlone = " << nlone << "; nbond = " << nbond << endl;
		}
		//assert(nlone == 2);
		//assert(nbond == 2);
	}//iwater
	return;
}

void Wan_centers_stat::calculate(const Water* water, const Cell &cel, int** &bond_wan_index, int** &lone_wan_index, const int &ito, const int &ith)
{
	for (int iwater=0; iwater<cel.atom[ito].na; iwater++)
	{
		double ang = HBs::angle(cel, cel.wan_centers[lone_wan_index[iwater][0]],
		cel.atom[ito].pos[iwater], cel.wan_centers[lone_wan_index[iwater][1]]);
		int which_angle = (int) (ang/INPUT.dtheta);
		this->lone_angle_distr[which_angle]++;
		// angle distribution done

		Vector3<double> pos_lone1 = cel.wan_centers[lone_wan_index[iwater][0]] - cel.atom[ito].pos[iwater];
		Vector3<double> pos_lone2 = cel.wan_centers[lone_wan_index[iwater][1]] - cel.atom[ito].pos[iwater];

		put_back_cell(pos_lone1);
		put_back_cell(pos_lone2);
		//cout << pos_lone1.x << " " << pos_lone1.y << " " << pos_lone1.z << endl;
		//cout << pos_lone2.x << " " << pos_lone2.y << " " << pos_lone2.z << endl;

		Vector3<double> pos_lone1_norm = pos_lone1/pos_lone1.norm();
		Vector3<double> pos_lone2_norm = pos_lone2/pos_lone2.norm();

		Vector3<double> yaxis = (pos_lone2_norm+pos_lone1_norm)/(pos_lone2_norm+pos_lone1_norm).norm();
		Vector3<double> xaxis = (pos_lone2_norm-pos_lone1_norm)/(pos_lone2_norm-pos_lone1_norm).norm();

		double x1 = pos_lone1.x*xaxis.x+
					pos_lone1.y*xaxis.y+
					pos_lone1.z*xaxis.z;
		
		double y1 = pos_lone1.x*yaxis.x+
					pos_lone1.y*yaxis.y+
					pos_lone1.z*yaxis.z;

		double x2 = pos_lone2.x*xaxis.x+
					pos_lone2.y*xaxis.y+
					pos_lone2.z*xaxis.z;
		
		double y2 = pos_lone2.x*yaxis.x+
					pos_lone2.y*yaxis.y+
					pos_lone2.z*yaxis.z;

		x1 += INPUT.rcut/2;
		x2 += INPUT.rcut/2;

		if (x1 < 0 or x2 < 0 or y1 < 0 or y2 < 0)
		{
			cout << "Warning! rcut is too small to cover Wannier centers" << endl;
			cout << "rcut = " << INPUT.rcut << "; x1 = " << x1 << "; x2 = " << x2 << endl;
			cout << "y1 = " << y1 << "; y2 = " << y2 << endl;
			exit(0);
		}
		if (x1 > INPUT.rcut or x2 > INPUT.rcut or y1 > INPUT.rcut or y2 > INPUT.rcut)
		{
			cout << "Warning! rcut is too small to cover Wannier centers" << endl;
			cout << "rcut = " << INPUT.rcut << "; x1 = " << x1 << "; x2 = " << x2 << endl;
			cout << "y1 = " << y1 << "; y2 = " << y2 << endl;
			exit(0);
		}
		int which_x1 = (int) (x1/INPUT.dr);
		int which_y1 = (int) (y1/INPUT.dr);
		int which_x2 = (int) (x2/INPUT.dr);
		int which_y2 = (int) (y2/INPUT.dr);
		this->lone_2D_distr[which_x1][which_y1]++;
		this->lone_2D_distr[which_x2][which_y2]++;
		
	}// iwater
	
	return;
}

void Wan_centers_stat::out(int &count_geometry_number)
{
	for (int ir=0; ir<this->nr; ir++)
	{
		for(int ir2=0; ir2<this->nr; ir2++)
		{
			this->lone_2D_distr[ir][ir2] /= count_geometry_number;
		}//ir2
	}//ir
	for(int iangle=0; iangle<this->nangle; iangle++)
	{
		this->lone_angle_distr[iangle] /= count_geometry_number;
	}// iangle

	double sum = 0;
	for (int ir=0; ir<this->nr; ir++)
	{
		for(int ir2=0; ir2<this->nr; ir2++)
		{
			sum += this->lone_2D_distr[ir][ir2]*INPUT.dr*INPUT.dr;
		}//ir2
	}//ir
	for (int ir=0; ir<this->nr; ir++)
	{
		for(int ir2=0; ir2<this->nr; ir2++)
		{
			this->lone_2D_distr[ir][ir2] /= sum;
		}//ir2
	}//ir

	sum = 0;
	for(int iangle=0; iangle<this->nangle; iangle++)
	{
		sum += this->lone_angle_distr[iangle]*INPUT.dtheta;
	}
	for(int iangle=0; iangle<this->nangle; iangle++)
	{
		this->lone_angle_distr[iangle] /= sum;
	}
	// end normalization
	ofstream ofs_angle("lone_wan_angle.txt");
	for(int iangle=0; iangle<this->nangle; iangle++)
	{
		ofs_angle << iangle*INPUT.dtheta+0.5*INPUT.dtheta << " " << this->lone_angle_distr[iangle] << endl;
	}
	ofs_angle.close();

	ofstream ofs_2D("lone_wan_2D_distribution.txt");
	ofs_2D << " ";
	for(int ix=0; ix<this->nr; ix++)
	{
		ofs_2D << " " << 0.5*INPUT.dr + ix*INPUT.dr;
	}
	ofs_2D << endl;
	for(int ix=0; ix<this->nr; ix++)
	{
		ofs_2D << 0.5*INPUT.dr + ix*INPUT.dr << " ";
		for(int iy=0; iy<this->nr; iy++)
		{
			ofs_2D << this->lone_2D_distr[ix][iy] << " ";
		}
		ofs_2D << endl;
	}
	ofs_2D.close();
	return;
}

void Wan_centers_stat::put_back_cell(Vector3<double> &pos)
{
	if (abs(pos.x) >= INPUT.celldm1/2)
	{
		if (abs(pos.x-INPUT.celldm1) < abs(pos.x) and abs(pos.x-INPUT.celldm1) < abs(pos.x+INPUT.celldm1))
		{
			pos.x -= INPUT.celldm1;
		}
		else if (abs(pos.x+INPUT.celldm1) < abs(pos.x) and abs(pos.x+INPUT.celldm1) < abs(pos.x-INPUT.celldm1))
		{
			pos.x += INPUT.celldm1;
		}
	}
	if (abs(pos.y) >= INPUT.celldm2/2)
	{
		if (abs(pos.y-INPUT.celldm2) < abs(pos.y) and abs(pos.y-INPUT.celldm2) < abs(pos.y+INPUT.celldm2))
		{
			pos.y -= INPUT.celldm2;
		}
		else if (abs(pos.y+INPUT.celldm2) < abs(pos.y) and abs(pos.y+INPUT.celldm2) < abs(pos.y-INPUT.celldm2))
		{
			pos.y += INPUT.celldm2;
		}
	}
	if (abs(pos.z) >= INPUT.celldm3/2)
	{
		if (abs(pos.z-INPUT.celldm3) < abs(pos.z) and abs(pos.z-INPUT.celldm3) < abs(pos.z+INPUT.celldm3))
		{
			pos.z -= INPUT.celldm3;
		}
		else if (abs(pos.z+INPUT.celldm3) < abs(pos.z) and abs(pos.z+INPUT.celldm3) < abs(pos.z-INPUT.celldm3))
		{
			pos.z += INPUT.celldm3;
		}
	}
	return;
}