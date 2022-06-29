#include "wannier1.h"
#include "input.h"
#include "cellFile.h"
#include "HBs.h"

Wannier1::Wannier1() 
{
}

Wannier1::~Wannier1() 
{
}

void Wannier1::Routine()
{
	// 2D data
    assert(INPUT.nx>0);
    assert(INPUT.ny>0);

    this->nx = INPUT.nx;
    this->ny = INPUT.ny;
	this->nz = INPUT.nz;
	this->x0 = INPUT.x0; // density [0.90,1.12] in g/mL
	this->y0 = INPUT.y0; // length  [0.90,1.06] in Angstroms
	this->z0 = INPUT.z0;
	this->dx = INPUT.dx; // 0.004
	this->dy = INPUT.dy; // 0.004
	this->dz = INPUT.dz;

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

	// mohan added 2018-08-20
	this->data_xyz = new double**[nx];
	for(int ix=0; ix<nx; ++ix)
	{
		this->data_xyz[ix] = new double*[ny];
		for(int iy=0; iy<ny; ++iy)
		{
			this->data_xyz[ix][iy] = new double[nz]();
		}
	}


	// plot out the final distance between H and wannier centers
	this->dr=INPUT.dr;
	assert(INPUT.dr>0);
	double rcut=INPUT.rcut;
	this->nr=rcut/dr+10; //+10 for safety
	dis_Hwan = new double[nr]();

	// setup geometry index
	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

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
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "snapshot " << cel.snapshot_index << " igeo " << igeo << endl;

		//read_wan(cel, dis_Hwan);
		//angle_analysis(cel, dis_Hwan);
		analysis2(cel);
		//analysis3(cel);
	}	

	// print out 1D Wannier centers
	ofstream ofs("dis_1D.dat");

	double sum=0.0;
	for(int ir=0; ir<nr; ++ir)
	{
		sum += dis_Hwan[ir]*dr;
	}
	if(sum>0.0)
	{
		for(int ir=0; ir<nr; ++ir)
		{
			dis_Hwan[ir]/=sum;
		}
	}
	sum=0.0;
	for(int ir=0; ir<nr; ++ir)
	{
		sum += dis_Hwan[ir]*dr;
	}
	cout << "sum is " << sum << endl;

	for(int ir=0; ir<nr; ++ir)
	{
		ofs << dr*(ir+1) << " " << dis_Hwan[ir] << endl;
	}
	ofs.close();
	delete[] dis_Hwan;


	// print out the information
	ofstream ofs2D("dis_2D.dat");
	
    sum=0.0;
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

	// print out the information
	ofstream ofs3D("dis_3D.dat");
	sum=0.0;
	for(int iz=0; iz<nz; ++iz)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			for(int ix=0; ix<nz; ++ix)
			{
				sum += data_xyz[ix][iy][iz] * INPUT.dx * INPUT.dy * INPUT.dz;
			}
		}
	}

	if(sum>0.0)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			for(int iy=0; iy<ny; ++iy)
			{
				for(int iz=0; iz<nz; ++iz)
				{
					//if(data_xyz[ix][iy][iz]>0)
					{
						ofs3D << ix*dx << " " << iy*dy << " " << iz*dz << " " << data_xyz[ix][iy][iz]/sum << endl;
					}
				}
			}
		}
	}
	ofs3D.close();

	// clean
    for(int ix=0; ix<nx; ++ix)
    {
        delete[] coord_xy[ix];
    }
    delete[] coord_xy;

	for(int ix=0; ix<nx; ++ix)
	{
		for(int iy=0; iy<ny; ++iy)
		{
			delete[] data_xyz[ix][iy];
		}
		delete[] data_xyz[ix];
	}
	delete[] data_xyz;

	return;
}

void Wannier1::analysis3(Cell &cel)
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
    if(INPUT.ntype==3)
    {
        cout << "ito = " << ito << endl;
        cout << "ith = " << ith << endl;
        cout << "itc = " << itc << endl;
        assert(itc>=0);
    }

    const double norm1 = cel.a1.norm();
    const double norm2 = cel.a2.norm();
    const double norm3 = cel.a3.norm();

    Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;
    HBs::setup_water(cel, water);

    int it3=ito;

    for(int ia=0; ia<cel.atom[ito].na; ++ia)
    {
        if(water[ia].nH!=2) continue;

        const int ih1 = water[ia].indexH[0];
        const int ih2 = water[ia].indexH[1];

		for(int ib=0; ib<INPUT.nbands; ++ib)
		{
			// distance between Wannier centers and oxygen atoms
			double dx1=shortest(cel.atom[ito].pos[ia].x, cel.wan_centers[ib].x, norm1);	
			double dy1=shortest(cel.atom[ito].pos[ia].y, cel.wan_centers[ib].y, norm2);	
			double dz1=shortest(cel.atom[ito].pos[ia].z, cel.wan_centers[ib].z, norm3);
			double dis1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);

			int indexZ = (dis1-z0)/dz;
	
			if(dis1>0.42) continue;

			for(int ia2=0; ia2<cel.atom[it3].na; ++ia2)
			{
				if(it3==ito and water[ia2].nH!=2) continue;
				if(it3==ito and ia==ia2) continue;

				for(int ia3=0; ia3<cel.atom[ith].na; ++ia3)
				{
					double dx2=shortest(cel.atom[it3].pos[ia2].x, cel.atom[ith].pos[ia3].x, norm1);
					double dy2=shortest(cel.atom[it3].pos[ia2].y, cel.atom[ith].pos[ia3].y, norm2);
					double dz2=shortest(cel.atom[it3].pos[ia2].z, cel.atom[ith].pos[ia3].z, norm3);
					double dis2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
					if(dis2>1.2) continue;

					double dx3=shortest(cel.atom[ito].pos[ia].x, cel.atom[ith].pos[ia3].x, norm1);
					double dy3=shortest(cel.atom[ito].pos[ia].y, cel.atom[ith].pos[ia3].y, norm2);
					double dz3=shortest(cel.atom[ito].pos[ia].z, cel.atom[ith].pos[ia3].z, norm3);
					double dis3 = sqrt(dx3*dx3+dy3*dy3+dz3*dz3);

					int indexX = (dis3-x0)/dx;

					double a3=angle(cel, cel.atom[ith].pos[ia3], cel.atom[ito].pos[ia], cel.wan_centers[ib]);

					int indexY = (a3-y0)/dy;

					if(indexX<nx and indexY<ny and indexZ<nz and indexX>=0 and indexY>=0 and indexZ>=0)
					{
						data_xyz[indexX][indexY][indexZ]+=1.0;
					}

				} // ia3

			} // ia2

		} // ib
	}


	delete[] water;
}


void Wannier1::analysis2(Cell &cel)
{
	// get ito, ith, and itc.
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
	if(INPUT.ntype==3)
	{ 
		cout << "ito = " << ito << endl;
		cout << "ith = " << ith << endl;
		cout << "itc = " << itc << endl;
		assert(itc>=0); 
	}

	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	// search for hexane	
	Water *water = new Water[cel.atom[ito].na];
	Water::nions = 0;
	HBs::setup_water(cel, water);

	int ion=-1;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		if(water[ia].nH==1) ion=ia;
	}
	cout << "ion is " << ion << endl;


	int it3=itc;

 	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		if(water[ia].nH!=2) continue;

		const int ih1 = water[ia].indexH[0]; 
		const int ih2 = water[ia].indexH[1];

			// mohan added
			double dx0=shortest(cel.atom[ito].pos[ia].x, cel.atom[ito].pos[ion].x, norm1);	
			double dy0=shortest(cel.atom[ito].pos[ia].y, cel.atom[ito].pos[ion].y, norm2);	
			double dz0=shortest(cel.atom[ito].pos[ia].z, cel.atom[ito].pos[ion].z, norm3);
			double dis0 = sqrt(dx0*dx0+dy0*dy0+dz0*dz0);
			if(dis0<12.0) continue;


		for(int ib=0; ib<INPUT.nbands; ++ib)
		{
			// distance between Wannier centers and oxygen atoms
			double dx1=shortest(cel.atom[ito].pos[ia].x, cel.wan_centers[ib].x, norm1);	
			double dy1=shortest(cel.atom[ito].pos[ia].y, cel.wan_centers[ib].y, norm2);	
			double dz1=shortest(cel.atom[ito].pos[ia].z, cel.wan_centers[ib].z, norm3);
			double dis1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
	
			if(dis1>0.42) continue;

			for(int ia2=0; ia2<cel.atom[it3].na; ++ia2)
			{
				if(it3==ito and water[ia2].nH!=2) continue;
				if(it3==ito and ia==ia2) continue;

				for(int ia3=0; ia3<cel.atom[ith].na; ++ia3)
				{
					double dx2=shortest(cel.atom[it3].pos[ia2].x, cel.atom[ith].pos[ia3].x, norm1);
					double dy2=shortest(cel.atom[it3].pos[ia2].y, cel.atom[ith].pos[ia3].y, norm2);
					double dz2=shortest(cel.atom[it3].pos[ia2].z, cel.atom[ith].pos[ia3].z, norm3);
					double dis2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
					if(dis2>1.2) continue;


					double dx3=shortest(cel.atom[ito].pos[ia].x, cel.atom[ith].pos[ia3].x, norm1);
					double dy3=shortest(cel.atom[ito].pos[ia].y, cel.atom[ith].pos[ia3].y, norm2);
					double dz3=shortest(cel.atom[ito].pos[ia].z, cel.atom[ith].pos[ia3].z, norm3);
					double dis3 = sqrt(dx3*dx3+dy3*dy3+dz3*dz3);

					// for C-O-H
					/*
					double dx2=shortest(cel.atom[ito].pos[ia].x, cel.atom[ith].pos[ia3].x, norm1);
					double dy2=shortest(cel.atom[ito].pos[ia].y, cel.atom[ith].pos[ia3].y, norm2);
					double dz2=shortest(cel.atom[ito].pos[ia].z, cel.atom[ith].pos[ia3].z, norm3);
					double dis2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
					if(dis2>1.2) continue;

					double dx3=shortest(cel.atom[ito].pos[ia].x, cel.atom[it3].pos[ia2].x, norm1);
					double dy3=shortest(cel.atom[ito].pos[ia].y, cel.atom[it3].pos[ia2].y, norm2);
					double dz3=shortest(cel.atom[ito].pos[ia].z, cel.atom[it3].pos[ia2].z, norm3);
					double dis3 = sqrt(dx3*dx3+dy3*dy3+dz3*dz3);
					*/

					int indexX = (dis3-x0)/dx;

/*					double a1=angle(cel, cel.atom[itc].pos[ia2], cel.atom[ito].pos[ia], cel.atom[ith].pos[ih1]);

					int indexY = (a1-y0)/dy;

					if(indexX<nx and indexY<ny and indexX>=0 and indexY>=0)
					{
						coord_xy[indexX][indexY]+=1.0;
					}

					double a2=angle(cel, cel.atom[itc].pos[ia2], cel.atom[ito].pos[ia], cel.atom[ith].pos[ih2]);
					indexY = (a2-y0)/dy;

					if(indexX<nx and indexY<ny and indexX>=0 and indexY>=0)
					{
						coord_xy[indexX][indexY]+=1.0;
					}
*/

					//double a3=angle(cel, cel.atom[ith].pos[ia3], cel.atom[ito].pos[ia], cel.wan_centers[ib]);
					//double a3=angle(cel, cel.atom[it3].pos[ia2], cel.atom[ito].pos[ia], cel.atom[ith].pos[ia3]);

					
					//int indexY = (a3-y0)/dy;
					int indexY = (dis1-y0)/dy;

					if(indexX<nx and indexY<ny and indexX>=0 and indexY>=0)
					{
						coord_xy[indexX][indexY]+=1.0;
					}

				} // ia3

			} // ia2

		} // ib
	}


	delete[] water;
}


void Wannier1::angle_analysis(Cell &cel, double* dis_Hwan)
{
	// get ito, ith, and itc.
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
	if(INPUT.ntype==3)
	{ 
		cout << "ito = " << ito << endl;
		cout << "ith = " << ith << endl;
		cout << "itc = " << itc << endl;
		assert(itc>=0); 
	}

	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	// search for hexane	
	Water *water = new Water[cel.atom[ito].na];
	Water::nions = 0;
	HBs::setup_water(cel, water);

 	for(int ia=0; ia<cel.atom[itc].na; ++ia)
	{
//		for(int ia4=0; ia4<cel.atom[ith].na; ++ia4)
//		{
//			double dx4=shortest(cel.atom[itc].pos[ia].x, cel.atom[ith].pos[ia4].x, norm1);
//			double dy4=shortest(cel.atom[itc].pos[ia].y, cel.atom[ith].pos[ia4].y, norm2);
//			double dz4=shortest(cel.atom[itc].pos[ia].z, cel.atom[ith].pos[ia4].z, norm3);
//			double dis4 = sqrt(dx4*dx4+dy4*dy4+dz4*dz4);

//			if(dis4<1.3)
//			{
				for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
				{
					if(water[ia2].nH!=2) continue;

					double dx2=shortest(cel.atom[itc].pos[ia].x, cel.atom[ito].pos[ia2].x, norm1);
					double dy2=shortest(cel.atom[itc].pos[ia].y, cel.atom[ito].pos[ia2].y, norm2);
					double dz2=shortest(cel.atom[itc].pos[ia].z, cel.atom[ito].pos[ia2].z, norm3);
					double dis2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);

					if(dis2<4)
					{
						   const int ih1 = water[ia2].indexH[0]; 
						   const int ih2 = water[ia2].indexH[1];

						   double a1=angle(cel, cel.atom[itc].pos[ia], cel.atom[ito].pos[ia2], cel.atom[ith].pos[ih1]);
						   double a2=angle(cel, cel.atom[itc].pos[ia], cel.atom[ito].pos[ia2], cel.atom[ith].pos[ih1]);
					//	   cout << "angle " << a1 << " " << a2 << endl;

						   int index=a1/this->dr;
						   dis_Hwan[index]++;

						   index=a2/this->dr;
						   dis_Hwan[index]++;

/*
						for(int ib=0; ib<INPUT.nbands; ++ib)
						{
							// distance between Wannier centers and oxygen atoms
							double dx3=shortest(cel.atom[ito].pos[ia2].x, cel.wan_centers[ib].x, norm1);	
							double dy3=shortest(cel.atom[ito].pos[ia2].y, cel.wan_centers[ib].y, norm2);	
							double dz3=shortest(cel.atom[ito].pos[ia2].z, cel.wan_centers[ib].z, norm3);
							double dis3 = sqrt(dx3*dx3+dy3*dy3+dz3*dz3);

							// for 1D
							if(dis3<0.42)
							{
								double a3=angle(cel, cel.atom[itc].pos[ia], cel.atom[ito].pos[ia2], cel.wan_centers[ib]);
								int index=a3/this->dr;
								dis_Hwan[index]++;
							}
						}
*/

					}
				}

//			} // 1.3
//		} // ith
	}// itc

	delete[] water;


	return;
}

double Wannier1::angle(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &pos3)
{
    double x1 = pos1.x; double y1 = pos1.y; double z1 = pos1.z;
    double x2 = pos2.x; double y2 = pos2.y; double z2 = pos2.z;
    double x3 = pos3.x; double y3 = pos3.y; double z3 = pos3.z;

    double d12x = shortest(x1, x2, INPUT.celldm1);
    double d12y = shortest(y1, y2, INPUT.celldm2);
    double d12z = shortest(z1, z2, INPUT.celldm3);

    double d32x = shortest(x3, x2, INPUT.celldm1);
    double d32y = shortest(y3, y2, INPUT.celldm2);
    double d32z = shortest(z3, z2, INPUT.celldm3);

    double norm1 = sqrt(d12x*d12x+d12y*d12y+d12z*d12z);
    double norm2 = sqrt(d32x*d32x+d32y*d32y+d32z*d32z);

    double angle = (d12x*d32x+d12y*d32y+d12z*d32z)/norm1/norm2;

    angle = acos(angle)/3.1415926535897*180;

    return angle;
}



void Wannier1::read_wan(Cell &cel, double* dis_Hwan)
{
	// get ito, ith, and itc.
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
	if(INPUT.ntype==3)
	{ 
		cout << "ito = " << ito << endl;
		cout << "ith = " << ith << endl;
		cout << "itc = " << itc << endl;
		assert(itc>=0); 
	}

	const double norm1 = cel.a1.norm();
	const double norm2 = cel.a2.norm();
	const double norm3 = cel.a3.norm();

	// search for hexane	
	Water *water = new Water[cel.atom[ito].na];
	Water::nions = 0;
	HBs::setup_water(cel, water);







	// hydroxide
	int ia_ion=-1;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		//cout << "oxygen " << ia << endl;
		bool found=false;
		for(int ia2=0; ia2<cel.atom[itc].na; ++ia2)
		{
			double dx2=shortest(cel.atom[ito].pos[ia].x, cel.atom[itc].pos[ia2].x, norm1);
			double dy2=shortest(cel.atom[ito].pos[ia].y, cel.atom[itc].pos[ia2].y, norm2);
			double dz2=shortest(cel.atom[ito].pos[ia].z, cel.atom[itc].pos[ia2].z, norm3);
			double dis2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);

			if(dis2<4.0)
			{
				found=true;
				for(int ib=0; ib<INPUT.nbands; ++ib)
				{
					double dx=shortest(cel.atom[ito].pos[ia].x, cel.wan_centers[ib].x, norm1);
					double dy=shortest(cel.atom[ito].pos[ia].y, cel.wan_centers[ib].y, norm2);
					double dz=shortest(cel.atom[ito].pos[ia].z, cel.wan_centers[ib].z, norm3);
					double dis = sqrt(dx*dx+dy*dy+dz*dz);
					if(dis<0.8)
					{
//						if(dis<0.31) 
//						{
							//cout << "water information " << water[ia].naccept << " " << water[ia].ndonate << endl;
//						}

						int index=dis/this->dr;
						dis_Hwan[index]++;
		//				cout << ia2 << " " << ib << endl;
					}
				}
			}
			if(found)break;
		}


		if(water[ia].nH==1)
		{
			ia_ion=ia;
			cout << "ia " << ia << " is ion" << endl;

			for(int ib=0; ib<INPUT.nbands; ++ib)
			{
				// distance between Wannier centers and oxygen atoms
				double dx=shortest(cel.atom[ito].pos[ia].x, cel.wan_centers[ib].x, norm1);	
				double dy=shortest(cel.atom[ito].pos[ia].y, cel.wan_centers[ib].y, norm2);	
				double dz=shortest(cel.atom[ito].pos[ia].z, cel.wan_centers[ib].z, norm3);
				double dis = sqrt(dx*dx+dy*dy+dz*dz);

				// for 1D
				if(dis<0.8) 
				{
					int index=dis/this->dr;
//					assert(index<this->nr);
	//				dis_Hwan[index]++;
				}
			}
		}
	}
	assert(ia_ion!=-1);



	// carbon and surronding water
	double dis_min0=100;
	int ia_o=-1;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
			const double dx1 = shortest(cel.atom[itc].pos[24].x, cel.atom[ito].pos[ia].x, norm1);
			const double dy1 = shortest(cel.atom[itc].pos[24].y, cel.atom[ito].pos[ia].y, norm2);
			const double dz1 = shortest(cel.atom[itc].pos[24].z, cel.atom[ito].pos[ia].z, norm3);
			double dis = sqrt(dx1*dx1+dy1*dy1+dz1*dz1); 

			const double dx2 = shortest(cel.atom[itc].pos[24].x, cel.atom[ito].pos[ia_ion].x, norm1);
			const double dy2 = shortest(cel.atom[itc].pos[24].y, cel.atom[ito].pos[ia_ion].y, norm2);
			const double dz2 = shortest(cel.atom[itc].pos[24].z, cel.atom[ito].pos[ia_ion].z, norm3);
			double dis2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2); 
		
			if(dis+dis2<dis_min0)
			{	
				dis_min0=dis;
				ia_o=ia;
			}
	}
	cout << "oxygen " << ia_o << " dis " << dis_min0 << endl;










	double dis_min=100;
	int ia_c=-1;
	for(int ia=0; ia<cel.atom[itc].na; ++ia)
	{
//		cout << "carbon " << ia << endl;
		for(int ia2=0; ia2<cel.atom[ith].na; ++ia2)
		{
			const double dx1 = shortest(cel.atom[itc].pos[ia].x, cel.atom[ith].pos[ia2].x, norm1);
			const double dy1 = shortest(cel.atom[itc].pos[ia].y, cel.atom[ith].pos[ia2].y, norm2);
			const double dz1 = shortest(cel.atom[itc].pos[ia].z, cel.atom[ith].pos[ia2].z, norm3);
			double dis = sqrt(dx1*dx1+dy1*dy1+dz1*dz1); 
			if(dis<1.3)
			{
				int index=dis/this->dr;
//				assert(index<this->nr);
//				dis_Hwan[index]++;
	//			cout << ia2+1 << " " << dis << endl;
			}
		}

		double dx2 = shortest(cel.atom[itc].pos[ia].x, cel.atom[ito].pos[ia_ion].x, norm1);
		double dy2 = shortest(cel.atom[itc].pos[ia].y, cel.atom[ito].pos[ia_ion].y, norm2);
		double dz2 = shortest(cel.atom[itc].pos[ia].z, cel.atom[ito].pos[ia_ion].z, norm3);
		double dis2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
		if(dis2 < dis_min)
		{
			dis_min = dis2;
			ia_c = ia;
		}

/*
		for(int ia3=0; ia3<cel.atom[itc].na; ++ia3)
		{
			const double dx1 = shortest(cel.atom[itc].pos[ia].x, cel.atom[itc].pos[ia3].x, norm1);
			const double dy1 = shortest(cel.atom[itc].pos[ia].y, cel.atom[itc].pos[ia3].y, norm2);
			const double dz1 = shortest(cel.atom[itc].pos[ia].z, cel.atom[itc].pos[ia3].z, norm3);
			double dis = sqrt(dx1*dx1+dy1*dy1+dz1*dz1); 
			if(dis<3.0)
			{
				cout << ia3 << " " << dis << endl;
			}

		}
		int ok;
		cin >> ok;
*/
	}
	cout << "carbon " << ia_c << " to O " << ia_ion << " distance " << dis_min << endl;








	//int it0=ito;
	int it0=itc;
	for(int ia=0; ia<cel.atom[it0].na; ++ia)
	{
        //if(ia==1)
	//	if(ia%6==5) // the xxx hexane
//		if(ia==24)
		{
//			cout << ia << endl;
			for(int ib=0; ib<INPUT.nbands; ++ib)
			{
				// distance between Wannier centers and oxygen atoms
				double dx=shortest(cel.atom[it0].pos[ia].x, cel.wan_centers[ib].x, norm1);	
				double dy=shortest(cel.atom[it0].pos[ia].y, cel.wan_centers[ib].y, norm2);	
				double dz=shortest(cel.atom[it0].pos[ia].z, cel.wan_centers[ib].z, norm3);
				double dis = sqrt(dx*dx+dy*dy+dz*dz);

				if(dis < 0.8) // 0.8 Angstroms is a safe cutoff
				{
					//cout << "C " << itc << " " << ia << " " << dis << endl;
				}

				// for 1D
				if(dis<INPUT.rcut) 
				{
					int index=dis/this->dr;
	//				assert(index<this->nr);
				//	dis_Hwan[index]++;
				}
			}// ib
		}
	}// itc


	delete[] water;

	return;
}
