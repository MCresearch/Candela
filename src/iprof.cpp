#include "cellFile.h"
#include "input.h"
#include "iprof.h"
#include "math.h"

void Iprofile::Routine()
{
	TITLE("Iprof","Routine");

	int nr = INPUT.iprof_nr;
	assert(nr>0);
	this->dr = 40.0/(double)nr;
	double b = INPUT.iprof_b;
	assert(b>0.0);
	double a = 1.0/sqrt(PI*2.0)/b;

	// unit is made for ionic density profile
	double* rab = new double[nr];
	for(int ir=0; ir<nr; ++ir)
	{
		rab[ir] = dr;
	}

	double* test = new double[nr]();
	for(int ir=0; ir<nr; ++ir)
	{
		double z = ir*dr-18;
		test[ir] = exp( -z*z/(b*b*2.0) ) * a; 
	}
	double tmp_sum=0.0;

	int tmp_nr = nr;
	if(nr%2==0) 
	{
      tmp_nr = nr-1;	
	}
    Math::Simpson_Integral(tmp_nr, test, rab, tmp_sum);
	cout << "tmp_sum= " << tmp_sum << endl;
	delete[] test;


	// allocate iprof
	double* r = new double[nr];
	double** iprof = new double*[INPUT.ntype];
	double** iprof_unit = new double*[INPUT.ntype];
	for(int it=0; it<INPUT.ntype; ++it)
	{
		iprof[it] = new double[nr];
		iprof_unit[it] = new double[nr];
		for(int ir=0; ir<nr; ++ir)
		{
			iprof[it][ir] = 0.0;
			iprof_unit[it][ir] = 0.0;
		}
	}


	assert(INPUT.geo_interval>0);
	int count_geometry_number=0;

	cout << "geo1=" << INPUT.geo_1 << endl;
	cout << "geo2=" << INPUT.geo_2 << endl;


	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;

		// cel : input geometry file
		
		CellFile cel;

		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;


		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;
		cout << "igeo=" << igeo << endl;

		// the original method
		//	cal( cel_in );

		// new method to calcultae the ionic density profile
		// 2014-03-09
		if(cel.read_and_used==false) continue;
		
		cal_gauss( cel, r, iprof, iprof_unit );
		++count_geometry_number;
	}


	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ir=0; ir<nr; ++ir)
		{
			iprof_unit[it][ir] = iprof[it][ir]/count_geometry_number;
		}
	}

	double sum = 0.0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ir=0; ir<nr; ++ir)
		{
			sum+=iprof_unit[it][ir]*dr;
		}
	}
	assert( sum > 0.0 );
	cout << " sum from ionic density profile = " << sum << endl;



	// print ionic density profile data
	
	if(RANK==0)
	{
		cout << " Print data to " << INPUT.geo_out << endl;
		ofstream ofs(INPUT.geo_out.c_str());
		ofs << "r_Ang ";
		
		for(int it=0; it<INPUT.ntype; ++it)
		{
			ofs << "type" << it+1 << " ";
		}
		ofs << endl;

		for(int ir=0; ir<nr; ++ir)
		{
			ofs << r[ir];
			for(int it=0; it<INPUT.ntype; ++it)
			{
				//ofs << " " << iprof[it][ir] << " " << iprof_unit[it][ir];
				ofs << " " << iprof_unit[it][ir]/surface_area;
			}
			ofs	<< endl; 
		}
		ofs.close();
	}

	// finish
	for(int it=0; it<INPUT.ntype; ++it)
	{
		delete[] iprof[it];
		delete[] iprof_unit[it];
	}
	delete[] iprof;
	delete[] iprof_unit;
	delete[] r;
	delete[] rab;

	return;
}

void Iprofile::cal_gauss( const Cell &cel, double* r, double** iprof, double** iprof_unit)
{
	TITLE("Iprofile","cal_gauss");

	// natom / v_liquid
	double norm1 = cel.a1.norm();
	double norm2 = cel.a2.norm();
	double norm3 = cel.a3.norm();
	assert(INPUT.iprof_nr>0);
	assert(norm1>0);
	assert(norm2>0);
	assert(norm3>0);

	this->surface_area = norm1 * norm2;

	// assume 'a3' is the longest axis
	this->dr = norm3/(double)INPUT.iprof_nr;
	for(int ir=0; ir<INPUT.iprof_nr; ++ir)
	{
		r[ir] = ir*dr;
	}

	// search for each atom
	assert( INPUT.iprof_b > 0.0 );
	double b = INPUT.iprof_b; // coefficient
	//cout << " coefficeint for gaussian broadening is " << b << endl; 

	double b2 = b * b;
	double a = 1.0/sqrt(2.0*PI)/b;
	int iat=0;

// only for carbon
	for(int it=0; it<INPUT.ntype; ++it)
	{
//		if(cel.atom[it].id!="O") continue;
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
//			if(INPUT.hindex==1){ if(ia!=0 and ia!=6  and ia!=12 and ia!=18 and ia!=24 and ia!=30) continue;}
//			else if(INPUT.hindex==2){ if(ia!=1 and ia!=7  and ia!=13 and ia!=19 and ia!=25 and ia!=31) continue;}
//			else if(INPUT.hindex==3){ if(ia!=2 and ia!=8  and ia!=14 and ia!=20 and ia!=26 and ia!=32) continue;}
//			else if(INPUT.hindex==4){ if(ia!=3 and ia!=9  and ia!=15 and ia!=21 and ia!=27 and ia!=33) continue;}
//			else if(INPUT.hindex==5){ if(ia!=4 and ia!=10 and ia!=16 and ia!=22 and ia!=28 and ia!=34) continue;}
//			else if(INPUT.hindex==6){ if(ia!=5 and ia!=11 and ia!=17 and ia!=23 and ia!=29 and ia!=35) continue;}

/*
			if(INPUT.hindex==1){ if(ia!=0) continue;}
			else if(INPUT.hindex==2){ if(ia!=6) continue;}
			else if(INPUT.hindex==3){ if(ia!=12) continue;}
			else if(INPUT.hindex==4){ if(ia!=18) continue;}
			else if(INPUT.hindex==5){ if(ia!=24) continue;}
			else if(INPUT.hindex==6){ if(ia!=30) continue;}
*/
			
			double posz = cel.atom[it].pos[ia].z;
			double dz = 0.0;


			while( posz >= norm3 )
			{
				posz -= norm3;	
			}
			while( posz < 0)
			{
				posz += norm3;
			}
			// calculate the gaussian distribution

			for(int ir=0; ir<INPUT.iprof_nr; ++ir)
			{
				double dz = posz - r[ir];
				for(int i=-1; i<=1; ++i)
				{
					double dz2 = dz + i*norm3;
					iprof[it][ir] += exp( -dz2*dz2/2.0/b2 )*a;
				}
			}
		}
	}


	return;
	// for oxygen
	for(int ia=0; ia<cel.atom[0].na; ++ia)
	{	
		int count=0;
		int indexH=-1;
		for(int ia2=0; ia2<cel.atom[2].na; ++ia2)
		{
			double dx = shortest(cel.atom[0].pos[ia].x, cel.atom[2].pos[ia2].x, norm1);
			double dy = shortest(cel.atom[0].pos[ia].y, cel.atom[2].pos[ia2].y, norm2);
			double dz = shortest(cel.atom[0].pos[ia].z, cel.atom[2].pos[ia2].z, norm3);
			double dis = sqrt(dx*dx+dy*dy+dz*dz);
			if(dis<1.24)
			{
				count++;
				indexH=ia2;
			}


		}
		//if(count==1) 
		if(count!=1)  //get oxygen(water) except hydroxide
		{
			double posz = cel.atom[0].pos[ia].z;
			while( posz >= norm3 )
			{
				posz -= norm3;	
			}
			while( posz < 0)
			{
				posz += norm3;
			}
			// calculate the gaussian distribution
			for(int ir=0; ir<INPUT.iprof_nr; ++ir)
			{
				double dz = posz - r[ir];
				double dz2 = dz * dz;
				double Gauss = exp( -dz2/b2 ) * a; 
				iprof[0][ir] += Gauss;
			}


			assert(indexH>=0);
			posz = cel.atom[2].pos[indexH].z;
			while( posz >= norm3 )
			{
				posz -= norm3;	
			}
			while( posz < 0)
			{
				posz += norm3;
			}
			// calculate the gaussian distribution
			for(int ir=0; ir<INPUT.iprof_nr; ++ir)
			{
				double dz = posz - r[ir];
				double dz2 = dz * dz;
				double Gauss = exp( -dz2/b2 ) * a; 
				iprof[2][ir] += Gauss;
			}
		}
	}

	return;
}
