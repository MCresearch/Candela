#include "density2D.h"
#include "input.h"
#include "cellFile.h"
#include "HBs.h"

Density2D::Density2D() 
{
}

Density2D::~Density2D() 
{
}

void Density2D::Routine()
{

    assert(INPUT.nx>0);
    assert(INPUT.ny>0);

    this->nx = INPUT.nx;
    this->ny = INPUT.ny;
	this->x0 = INPUT.x0; // density [0.90,1.12] in g/mL
	this->y0 = INPUT.y0; // length  [0.90,1.06] in Angstroms
	this->dx = INPUT.dx; // 0.004
	this->dy = INPUT.dy; // 0.004

	ofs_running << "nx = " << nx << endl;
	ofs_running << "ny = " << ny << endl;
	ofs_running << "x0 = " << x0 << endl;
	ofs_running << "y0 = " << y0 << endl;
	ofs_running << "dx = " << dx << endl;
	ofs_running << "dy = " << dy << endl;

	this->coord_xy = new double*[nx];
    for(int ix=0; ix<nx; ++ix)
    {
        this->coord_xy[ix] = new double[ny]();
    }


	// search for covalent bonds
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

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
		cout << "snapshot " << cel.snapshot_index << " igeo " << igeo << endl;
	
		covalent(cel, igeo);	
	}	

	// print out the information
	ofstream ofs2D("dis_cov_2D.dat");
	
    double sum=0.0;
    for(int iy=0; iy<ny; ++iy)
    {
        for(int ix=0; ix<nx; ++ix)
        {
            sum += coord_xy[ix][iy] * INPUT.dx * INPUT.dy;
        }
    }

    if(sum>0.0)
    {
        for(int iy=0; iy<ny; ++iy)
        {
            for(int ix=0; ix<nx; ++ix)
            {
                ofs2D << coord_xy[ix][iy]/sum << " ";
            }
            ofs2D << endl;
        }
    }

    for(int ix=0; ix<nx; ++ix)
    {
        delete[] coord_xy[ix];
    }
    delete[] coord_xy;

    ofs2D.close();
}

void Density2D::covalent(const Cell &cel, const int &igeo)
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

	double density = cel.atom[ito].na*18*1.6605/cel.volume;
	cout << "density=" << density << endl;

	int indexX = (density-x0)/dx;

	for(int o1=0; o1<cel.atom[ito].na; ++o1)
	{
		if(INPUT.func==1) // density-covalent relationship
		{
			for(int ih=0; ih<2; ++ih)
			{
				double dis=water[o1].disH[ih];
				int indexY = (dis-y0)/dy;
				if(indexX<nx and indexY<ny and indexX>=0 and indexY>=0)
				{
					coord_xy[indexX][indexY]+=1.0;
				}
			}
		}
		else if(INPUT.func==2) // H-bonded OO - covalent relationship
		{
			for(int io=0; io<water[o1].ndonate; ++io)
			{
				int indexO = water[o1].donateO[io];
				int indexH = water[o1].donateH[io];
			    double dis_x = distance(cel.atom[ito].pos[o1], cel.atom[ito].pos[indexO], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				double dis_y = distance(cel.atom[ito].pos[o1], cel.atom[ith].pos[indexH], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				int indexX2 = (dis_x-x0)/dx; 
				int indexY = (dis_y-y0)/dy;
				if(indexX2<nx and indexY<ny and indexX2>=0 and indexY>=0)
				{
					coord_xy[indexX2][indexY]+=1.0;
				}
			}
				
		}
		else if(INPUT.func==3) // H-bonded OO -- delta H
		{
			for(int io=0; io<water[o1].ndonate; ++io)
			{
				int indexO = water[o1].donateO[io];
				int indexH = water[o1].donateH[io];
			    double dis_x = distance(cel.atom[ito].pos[o1], cel.atom[ito].pos[indexO], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				double dis_y = distance(cel.atom[ito].pos[indexO], cel.atom[ith].pos[indexH], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3)
						-  distance(cel.atom[ito].pos[o1], cel.atom[ith].pos[indexH], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				int indexX2 = (dis_x-x0)/dx; 
				int indexY = (dis_y-y0)/dy;
				if(indexX2<nx and indexY<ny and indexX2>=0 and indexY>=0)
				{
					coord_xy[indexX2][indexY]+=1.0;
				}
			}
		}
	//	else if(INPUT.func==4) // delta H -- cos_theta
	}

	delete[] water;
	return;
}
