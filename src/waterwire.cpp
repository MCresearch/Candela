#include "cellFile.h"
#include "input.h"
#include "math.h"
#include "mj.h"
#include "HBs.h"
#include "waterwire.h"

Waterwire::Waterwire()
{
}

Waterwire::~Waterwire(){}

void Waterwire::Routine()
{
	TITLE("Waterwire","Routine");
	
	cout << "Compute the Properties of Waterwire" << endl;

	// setup bonds
	double dv = INPUT.bdf_dtheta;
	assert( dv > 0.0 );
	int npoints = 180.0/dv + 1;
	cout << " delta value of cos is = " << dv;
	cout << " number of points between -1 and 1 : " << npoints << endl;
	bond_df = new double[npoints](); // bond angle distribution function
	bond_df1 = new double[npoints]();
	bond_df2 = new double[npoints]();
	this->count_bonds=0;
	this->count_bonds1=0;
	this->count_bonds2=0;
	this->avg_bonds=0.0;
	this->avg_bonds1=0.0;
	this->avg_bonds2=0.0;
	
	// setup distances
	int nxx = INPUT.rcut/INPUT.dr;
	dis_o1 = new double[nxx]();
	dis_o2 = new double[nxx]();
	dis_o3 = new double[nxx]();
	dis_o123 = new double[nxx]();
	this->count_dis=0;
	this->avg_dis1=0.0;
	this->avg_dis2=0.0;
	this->avg_dis12=0.0;	

	// setup time-dependent variables
	this->time_nnn = 9600; // 600 for 0.5 ps
	this->time_dis1 = new double[time_nnn]();	
	this->time_dis1_max = new double[time_nnn]();
	this->time_dis1_min = new double[time_nnn]();
	this->time_count1 = new double[time_nnn]();	
	this->time_dis2 = new double[time_nnn]();
	this->time_dis2_max = new double[time_nnn]();
	this->time_dis2_min = new double[time_nnn]();
	this->time_count2 = new double[time_nnn]();	
	this->time_accept = new double[time_nnn]();
	for(int i=0; i<time_nnn; ++i)
	{
		time_dis1_min[i]=10000.0;
		time_dis2_min[i]=10000.0;
	}

	this->total_frames=0;

	// setup statistics for PT
	// typical 100 points with dt=0.005 ps
	assert(INPUT.npoints>0); // number of points
	assert(INPUT.ps_dt>0.0); // delta t in ps	
	this->np_tot = INPUT.npoints*2+1;
	this->pt_time = new double[np_tot]();	
	this->coor_num = new double[np_tot](); 

	// OPEN FILES 
	PT.setup_PT();

	double agv=0.0;
	int count_agv=0;
	for(int ip=1; ip<PT.npt; ++ip)
	{
		if(PT.type_pt[ip-1]!="rattling")
		{
			agv+= PT.snapshot_time[ip]-PT.snapshot_time_pt[ip-1];
			++count_agv;
		}
	}
	if(agv>0)
	{
		cout << "average is " << agv/count_agv << " ps" << endl;
	}
	
//	exit(0);

	assert(INPUT.nx>0);
	assert(INPUT.ny>0);
	int nx = INPUT.nx;
	int ny = INPUT.ny;

	this->coord_xy = new double*[nx];
	for(int ix=0; ix<nx; ++ix)
	{
		this->coord_xy[ix] = new double[ny]();
	}

	// BEGIN CALCULATING DATA
	this->count_geometry_number=0;
	int ipt=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
		if(igeo%INPUT.geo_interval!=0) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

//		if(cel.snapshot_time > 1.0) continue; //mohan test

		if(cel.read_and_used==false) continue;
		++count_geometry_number;
		cout << "igeo=" << igeo << endl;
		search_compression(cel, igeo);
	}	


	// print out data
	ofstream ofs_result("waterwire.dat");


	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			if(coord_xy[ix][iy]>0)
			{
				coord_xy[ix][iy]=-std::log(coord_xy[ix][iy]*INPUT.factor); 
			}
		}
	}

// check
//	cout << std::log10(10) << endl;
//	cout << std::log(10) << endl;
	

	for(int iy=0; iy<ny; ++iy)
	{
		for(int ix=0; ix<nx; ++ix)
		{
			ofs_result << coord_xy[ix][iy] << " "; 
		}
		ofs_result << endl;
	}
	ofs_result.close();


	// print the time-dependent variables
	ofstream ofs_time("time_o1o2.dat");
	ofs_time << "i dis1 min max dis2 min max acc" << endl;

	for(int i=0; i<time_nnn; ++i)
	{
		if(time_count1[i]>0 and time_count2[i]>0)
		{
			double aa = time_dis1[i]/(double)time_count1[i];
			double bb = time_dis2[i]/(double)time_count2[i];
			double cc = time_accept[i]/(double)time_count1[i];
			ofs_time << i << " " << aa
			<< " " << time_dis1_min[i] << " " << time_dis1_max[i]
			<< " " << bb 
			<< " " << time_dis2_min[i] << " " << time_dis2_max[i]
			<< " " << cc
			<< endl;
		}
	}	

	delete[] time_dis1;
	delete[] time_dis2;
	delete[] time_dis1_max;
	delete[] time_dis1_min;
	delete[] time_dis2_max;
	delete[] time_dis2_min;
	delete[] time_count1;
	delete[] time_count2;
	delete[] time_accept;
	ofs_time.close();


	//---------------------------------------------
	// plot the o1 o2 distance distribution
	
	if(count_dis>0)
	{
		ofs_running << "average dis o1: " << this->avg_dis1/(double)count_dis << endl;
		ofs_running << "average dis o2: " << this->avg_dis2/(double)count_dis << endl;
		ofs_running << "average dis o1+o2: " << this->avg_dis12/(double)count_dis << endl;
	}

	ofstream ofs_o1o2o3("dis_o1o2o3.dat");
	double sum1=0.0;
	double sum2=0.0;
	double sum3=0.0;
	double sum4=0.0;
	for(int i=0; i<nxx; ++i)
	{
		sum1 += dis_o1[i]*INPUT.dr;
		sum2 += dis_o2[i]*INPUT.dr;
		sum3 += dis_o3[i]*INPUT.dr;
		sum4 += dis_o123[i]*INPUT.dr;
	}
	if(sum1>0.0 and sum2>0.0 and sum3>0.0 and sum4>0.0)
	{
		for(int i=0; i<nxx; ++i)
		{
			dis_o1[i]/=sum1;	
			dis_o2[i]/=sum2;	
			dis_o3[i]/=sum3;	
			dis_o123[i]/=sum4;	
		}
	}
	for(int i=0; i<nxx; ++i)
	{
		double free_energy = -std::log(dis_o123[i]*sum4*INPUT.factor); 
		ofs_o1o2o3 << i*INPUT.dr << " " 
			<< dis_o123[i] << " " << dis_o123[i]*sum4*INPUT.factor << " " << free_energy << endl;


//		double free_energy = -std::log(dis_o2[i]*sum2*INPUT.factor); 
//		ofs_o1o2o3 << i*INPUT.dr << " " 
//			<< dis_o2[i] << " " << dis_o2[i]*sum2*INPUT.factor << " " << free_energy << endl;


//		double free_energy = -std::log(dis_o3[i]*sum3*INPUT.factor); 
//		ofs_o1o2o3 << i*INPUT.dr << " " 
//			<< dis_o3[i] << " " << dis_o3[i]*sum3*INPUT.factor << " " << free_energy << endl;
	}
	delete[] dis_o1;
	delete[] dis_o2;
	delete[] dis_o3;
	delete[] dis_o123;
	ofs_o1o2o3.close();

	//---------------------------------------------
	// plot the angle distribution
	
	if(count_bonds>0)
	{
		ofs_running << "average bond angle ooo: " << avg_bonds/(double)count_bonds << endl;
	}
	
	if(count_bonds1>0)
	{
		ofs_running << "average bond angle h1oo: " << avg_bonds1/(double)count_bonds1 << endl;
	}

	if(count_bonds2>0)
	{
		ofs_running << "average bond angle h2oo: " << avg_bonds2/(double)count_bonds2 << endl;
	}

	
	ofstream ofs_angle("angle.dat");
	ofs_angle << "angle ooo h1oo h2oo hoo" << endl;

	double sum = 0.0;
	sum1 = 0.0;
	sum2 = 0.0;
	for(int i=0; i<npoints; ++i)
	{
		sum += bond_df[i]*dv;
		sum1 += bond_df1[i]*dv;
		sum2 += bond_df2[i]*dv;
	}

	if(sum>0.0)
	{
		for(int i=0; i<npoints; ++i)
		{
			bond_df[i]/=sum;
			bond_df1[i]/=sum1;
			bond_df2[i]/=sum2;
			ofs_angle << i*dv << " " << bond_df[i] << " " << bond_df1[i] 
			<< " " << bond_df2[i] << " " << 0.5*bond_df1[i]+0.5*bond_df2[i] << endl;
		}
	}

	delete[] bond_df;
	delete[] bond_df1;
	delete[] bond_df2;
	ofs_angle.close();

	// print out coor number change during pt
	ofstream ofs_cn("coord_number.dat");
	for(int i=0; i<np_tot; ++i)
	{
		if(coor_num[i]==0.0) pt_time[i]=1.0;
		ofs_cn << (i-INPUT.npoints)*INPUT.ps_dt << " " << coor_num[i]/pt_time[i] << endl; 
	}
	ofs_cn.close();

	// clean up
	for(int ix=0; ix<nx; ++ix)
	{
		delete[] coord_xy[ix];
	}
	delete[] coord_xy;

	// for coor number change during pt
	delete[] pt_time;
	delete[] coor_num;

	return;
}


void Waterwire::search_compression(const Cell &cel, const int &igeo)
{
	// get ito, ith, and itc.
	int ito=-1;
	int ith=-1;
	int itc=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
//		cout << cel.atom[it].id << endl;
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
	}
//	cout << "ito=" << ito << endl;
//	cout << "ith=" << ith << endl;
	if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
	if(INPUT.ntype==3){ assert(itc>=0); }


	Water *water = new Water[cel.atom[ito].na];
	Water::nions=0;

	HBs::setup_water(cel, water);
	ofs_running << "snapshot " << cel.snapshot_index << " nions = " << Water::nions << endl;

	if(INPUT.func<=1)
	{
		single_jump(cel, water, ito, ith);
	}
	else if(INPUT.func>=2)
	{
		double_jump(cel, water, ito, ith);
	}

	delete[] water;
}

void Waterwire::double_jump(const Cell &cel, Water *water, const int &ito, const int &ith)
{
	total_frames=0;
	for(int ip=1; ip<PT.npt-2; ++ip)
	{
//		cout << "ip=" << ip << endl;
		if(PT.type_pt[ip-1]=="double")
		{
			++total_frames;

			// PT.snapshot_index_pt[ip-1]: time for h1 (first PT) to leave o1 
			// PT.snapshot_index[ip]: time for h1 to arrive o2
			// PT.snapshot_index_pt[ip]: time for h2 to leave o2
			// PT.snapshot_index[ip+1]: time for h2 to arrive o3
			int s1=PT.snapshot_index_pt[ip-1];
			int s2=PT.snapshot_index[ip+1];

			int dtt= s2-s1; 
			int dsmin = s1-dtt*INPUT.ext_1-INPUT.movement_x;
			int dsmax = s2+dtt*INPUT.ext_2+INPUT.movement_y;


//			if(cel.snapshot_index > dsmin and cel.snapshot_index < dsmax)
			{

				// indexes for atoms on the waterwire
				int o1 = PT.ion_index[ip-1]-1;
				int o2 = PT.ion_index[ip]-1;
				int o3 = PT.ion_index[ip+1]-1;
				int h1 = PT.indexH[ip]-1;
				int h2 = PT.indexH[ip+1]-1;
	

			// mohan added 2017-03-25
			/*
			if(INPUT.system=="hydronium")
			{
				if(water[o1].nH!=3 and water[o2].nH!=3 and water[o3].nH!=3)
				{
					continue;
				}
			}
*/

				// distance between o1 and o2, and o2 and o3
				double dis1 = distance(cel.atom[ito].pos[o1], cel.atom[ito].pos[o2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
				double dis4 = distance(cel.atom[ito].pos[o2], cel.atom[ito].pos[o3], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

				if(h1==-1) continue;
				if(h2==-1) continue;

				assert(o1>=0);
				assert(o2>=0);
				assert(o3>=0);
				assert(h1>=0);
				assert(h2>=0);
				
				// distances between o1 and h1, o2 and h1
				double dis2 = distance(cel.atom[ito].pos[o1], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
				double dis3 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

				// distances between o2 and h2, o3 and h2
				double dis5 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
				double dis6 = distance(cel.atom[ito].pos[o3], cel.atom[ith].pos[h2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

				// distance between o1 and o3
				double dis7 = distance(cel.atom[ito].pos[o1], cel.atom[ito].pos[o3], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 

				// angle between o1-o2-o3
				double angle = HBs::angle(cel, cel.atom[ito].pos[o1], cel.atom[ito].pos[o2], cel.atom[ito].pos[o3]);

				// angle formed by h1-o1-o2
				double angle1 = HBs::angle(cel, cel.atom[ith].pos[h1], cel.atom[ito].pos[o1], cel.atom[ito].pos[o2]);

				// angle formed by h2-o2-o3
				double angle2 = HBs::angle(cel, cel.atom[ith].pos[h2], cel.atom[ito].pos[o2], cel.atom[ito].pos[o3]);

				// waterwire length
				double r123 = dis1+dis4;

				// proton transfer coordinate 1 and 2
				double v1 = dis2 - dis3;
				double v2 = dis5 - dis6;
				double ptcoord = (v1 + v2)/2.0;

				// extremely important !!! mohan added 2017-01-18
/*
				bool should_count=false;
				if(cel.snapshot_index < s1 and v1 > 0 and v2 > 0) should_count=true;
				if(cel.snapshot_index < s1 and v1 < 0 and v2 < 0) should_count=true;
				if(cel.snapshot_index > s2 and v1 > 0 and v2 > 0) should_count=true;
				if(cel.snapshot_index > s2 and v1 < 0 and v2 < 0) should_count=true;
				if(cel.snapshot_index >= s1 and cel.snapshot_index<=s2) should_count=true;
*/

				bool should_count=true;
				if(INPUT.system=="hydronium")
				{
					// before first jump
					double tol=0.0;
					if(cel.snapshot_index < s1)
					{
						if(v1-tol>0 or v2-tol>0) should_count=false;
//						if(dis2>INPUT.rcut_oh) should_count=false;
//						if(dis5>INPUT.rcut_oh) should_count=false;
						//if(water[o1].nH!=3) should_count=false;
					}
					// after second jump
					if(cel.snapshot_index > s2)
					{
						if(v1+tol<0 or v2+tol<0) should_count=false;
//						if(dis3>INPUT.rcut_oh) should_count=false;
//						if(dis6>INPUT.rcut_oh) should_count=false;
						//if(water[o3].nH!=3) should_count=false;
					}

					// mohan added 2017-03-24
					bool formed_HB12 = false;
					bool formed_HB21 = false;
					bool formed_HB23 = false;
					bool formed_HB32 = false;
					for(int id=0; id<water[o1].ndonate; ++id) 
					if( water[o1].donateO[id] == o2 and water[o1].donateH[id] == h1) formed_HB12=true;
					for(int id=0; id<water[o2].ndonate; ++id) 
					if( water[o2].donateO[id] == o1 and water[o2].donateH[id] == h1) formed_HB21=true;
					for(int id=0; id<water[o2].ndonate; ++id) 
					if( water[o2].donateO[id] == o3 and water[o2].donateH[id] == h2) formed_HB23=true;
					for(int id=0; id<water[o3].ndonate; ++id) 
					if( water[o3].donateO[id] == o2 and water[o3].donateH[id] == h2) formed_HB32=true;
					if(formed_HB12 and formed_HB23){} 
					else if(formed_HB12 and formed_HB32){}
					else if(formed_HB21 and formed_HB32){}
					else if(formed_HB21 and formed_HB23){}
					else should_count=false;
				}
				if(INPUT.system=="hydroxide")
				{
					if(cel.snapshot_index < s1 and (v1 < 0 or v2 < 0))
					{
						should_count=false;
					}
					if(cel.snapshot_index > s2 and (v1 > 0 or v2 > 0))
					{
						should_count=false;
					}

					// mohan added 2017-03-24
					bool formed_HB12 = false;
					bool formed_HB21 = false;
					bool formed_HB23 = false;
					bool formed_HB32 = false;
					for(int id=0; id<water[o1].naccept; ++id) 
					if( water[o1].acceptO[id] == o2 and water[o1].acceptH[id] == h1) formed_HB12=true;
					for(int id=0; id<water[o2].naccept; ++id) 
					if( water[o2].acceptO[id] == o1 and water[o2].acceptH[id] == h1) formed_HB21=true;
					for(int id=0; id<water[o2].naccept; ++id) 
					if( water[o2].acceptO[id] == o3 and water[o2].acceptH[id] == h2) formed_HB23=true;
					for(int id=0; id<water[o3].naccept; ++id) 
					if( water[o3].acceptO[id] == o2 and water[o3].acceptH[id] == h2) formed_HB32=true;
					if(formed_HB12 and formed_HB23){} 
					else if(formed_HB12 and formed_HB32){}
					else if(formed_HB21 and formed_HB32){}
					else if(formed_HB21 and formed_HB23){}
					else should_count=false;
				}

				int indexX=-1;
				int indexY=-1;

				// typical 'arm' figure.
				if(INPUT.func==2)
				{
//					should_count = true; // just for test; renxi 20200507
					if(INPUT.system=="hydronium")
					{
						indexX = (ptcoord-INPUT.x0)/INPUT.dx;
					}
					if(INPUT.system=="hydroxide")
					{
						indexX = (-ptcoord-INPUT.x0)/INPUT.dx;
					}
					indexY = (r123-INPUT.y0)/INPUT.dy;
//					indexY = (dis7-INPUT.y0)/INPUT.dy;
				}
				// o1-o2 distance relation
				else if(INPUT.func==3)
				{
					indexX = (dis1-INPUT.x0)/INPUT.dx;
					indexY = (dis4-INPUT.y0)/INPUT.dy;
				}
				// ptcoordination - angle
				else if(INPUT.func==4)
				{
					indexX = (ptcoord-INPUT.x0)/INPUT.dx;
					indexY = (0.5*angle1+0.5*angle2-INPUT.y0)/INPUT.dy;
					//indexY = (angle1-INPUT.y0)/INPUT.dy;
				}
				else if(INPUT.func==5)
				{
					// x axis becomes the 'trajectory time'
					double tt1=PT.snapshot_time_pt[ip-1];
					double tt2=PT.snapshot_time[ip+1];
					indexX = (cel.snapshot_time - tt1 - INPUT.x0)/INPUT.dx;	
					// water wire length
					indexY = (r123-INPUT.y0)/INPUT.dy;
				}
				else if(INPUT.func==6)
				{
					if(INPUT.system=="hydronium")
					{
						indexX = (v1-INPUT.x0)/INPUT.dx;
					}
					indexY = (dis1-INPUT.y0)/INPUT.dy;
				}

			
				if(indexX<INPUT.nx and indexY<INPUT.ny and indexX>=0 and indexY>=0 and should_count==true)
				{
					int indexs = (cel.snapshot_index-dsmin)/10; // for PBE and PBE+vdW
					// int indexs = (cel.snapshot_index-dsmin)/5; // for PBE0+vdW
					if(indexs<time_nnn)
					{
						time_dis1[indexs] += dis1;
						if(dis1>time_dis1_max[indexs]) time_dis1_max[indexs]=dis1;
						if(dis1<time_dis1_min[indexs]) time_dis1_min[indexs]=dis1;
						time_count1[indexs]++;
						time_dis2[indexs] += dis4;
						if(dis4>time_dis2_max[indexs]) time_dis2_max[indexs]=dis4;
						if(dis4<time_dis2_min[indexs]) time_dis2_min[indexs]=dis4;
						time_count2[indexs]++;

						if(water[o1].naccept==1 and water[o1].ndonate==3) time_accept[indexs]+=1.0;
						else if(water[o2].naccept==1 and water[o2].ndonate==3) time_accept[indexs]+=1.0;
						else if(water[o3].naccept==1 and water[o3].ndonate==3) time_accept[indexs]+=1.0;
					}

					// bond angles
					int iii = int(angle/INPUT.bdf_dtheta); 
					this->bond_df[iii]+=1.0;
					this->avg_bonds+=angle;
					this->count_bonds++;

					if(dis2<INPUT.rcut_oh)
					{
						int iii1 = int(angle1/INPUT.bdf_dtheta); 
						this->bond_df1[iii1]+=1.0;
						this->avg_bonds1+=angle1;
						this->count_bonds1++;
					}

					if(dis5<INPUT.rcut_oh)
					{
						int iii2 = int(angle2/INPUT.bdf_dtheta); 
						this->bond_df2[iii2]+=1.0;
						this->avg_bonds2+=angle2;
						this->count_bonds2++;
					}

					// average distance O1-O2 and O2-O3
					if(dis1<INPUT.rcut and dis4<INPUT.rcut and dis7<INPUT.rcut and r123<INPUT.rcut)
					{
					//	cout << dis1 << " " << dis4 << " " << r123 << endl;
						int jjj = int(dis1/INPUT.dr);
						int kkk = int(dis4/INPUT.dr);
						int mmm = int(dis7/INPUT.dr);
						int lll = int(r123/INPUT.dr);
						this->dis_o1[jjj]+=1.0;
						this->dis_o2[kkk]+=1.0;
						this->dis_o3[mmm]+=1.0;
						this->dis_o123[lll]+=1.0;
						this->avg_dis1+=dis1;
						this->avg_dis2+=dis4;
						this->avg_dis12+=dis4+dis1;
						this->count_dis++;
					}
				
		
					if(INPUT.func_b==1) coord_xy[indexX][indexY]+=1.0;
					else if(INPUT.func_b==10)
					{
						if(water[o1].naccept==0) coord_xy[indexX][indexY]+=1.0;	
						else if(water[o2].naccept==0) coord_xy[indexX][indexY]+=1.0;	
						else if(water[o3].naccept==0) coord_xy[indexX][indexY]+=1.0;	
					}
					else if(INPUT.func_b==11)
					{
						if(water[o1].naccept==1) coord_xy[indexX][indexY]+=1.0;
						else if(water[o2].naccept==1) coord_xy[indexX][indexY]+=1.0;
						else if(water[o3].naccept==1) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==21)
					{
						if(water[o1].naccept==2) coord_xy[indexX][indexY]+=1.0;
						else if(water[o2].naccept==2) coord_xy[indexX][indexY]+=1.0;
						else if(water[o3].naccept==2) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==31)
					{
						if(water[o1].naccept==3) coord_xy[indexX][indexY]+=1.0;
						else if(water[o2].naccept==3) coord_xy[indexX][indexY]+=1.0;
						else if(water[o3].naccept==3) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==310)
					{
						if(water[o1].naccept==3 and water[o1].ndonate==0) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==311)
					{
						if(water[o1].naccept==3 and water[o1].ndonate==1) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==32)
					{
						if(water[o2].naccept==3) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==33)
					{
						if(water[o3].naccept==3) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==41)
					{
						if(water[o1].naccept>=4) coord_xy[indexX][indexY]+=1.0;
						else if(water[o2].naccept>=4) coord_xy[indexX][indexY]+=1.0;
						else if(water[o3].naccept>=4) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==410)
					{
						if(water[o1].naccept==4 and water[o1].ndonate==0) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==411)
					{
						if(water[o1].naccept==4 and water[o1].ndonate==1) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==42)
					{
						if(water[o2].naccept==4) coord_xy[indexX][indexY]+=1.0;
					}
					else if(INPUT.func_b==43)
					{
						if(water[o3].naccept==4) coord_xy[indexX][indexY]+=1.0;
					}

//					if(water[o1].naccept>=4 or water[o2].naccept>=4)
//					{
//						coord_xy[indexX][indexY]+=1.0;
//					}
//					if(water[o1].naccept==3 or water[o2].naccept==3)
//					{
//						coord_xy[indexX][indexY]+=1.0;
//					}
				}
			

			}
		}
	}
	return;
}


void Waterwire::single_jump(const Cell &cel, Water *water, const int &ito, const int &ith)
{
	//-------------------------------------------------
	// all pt events included, including rattling
	//-------------------------------------------------
	for(int ip=1; ip<PT.npt; ++ip)
	{
// please DIY
		if(PT.type_pt[ip-1]!="rattling")
		//if(PT.type_pt[ip-1]=="double")
		{
			// PT.snapshot_index_pt[ip-1]: time for h1 (first PT) to leave o1 
			// PT.snapshot_index[ip]: time for h1 to arrive o2
			// PT.snapshot_index_pt[ip]: time for h2 to leave o2
			// PT.snapshot_index[ip+1]: time for h2 to arrive o3
			int s1=PT.snapshot_index_pt[ip-1];
			int s2=PT.snapshot_index[ip];

			int o1 = PT.ion_index[ip-1]-1;
			int o2 = PT.ion_index[ip]-1;
			int h1 = PT.indexH[ip-1]-1;
			if(h1==-1) continue;

			assert(o1>=0);
			assert(o2>=0);
			assert(h1>=0);

			// distance between o1 and o2
			double dis1 = distance(cel.atom[ito].pos[o1], cel.atom[ito].pos[o2], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
			// distance between o1 and h1
			double dis2 = distance(cel.atom[ito].pos[o1], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
			// distance between o2 and h1
			double dis3 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h1], cel.a1.norm(), cel.a2.norm(), cel.a3.norm()); 
			// angle formed h1-o1-o2
			double angle1 = HBs::angle(cel, cel.atom[ith].pos[h1], cel.atom[ito].pos[o1], cel.atom[ito].pos[o2]);
			// angle formed h1-o2-o1
			double angle2 = HBs::angle(cel, cel.atom[ith].pos[h1], cel.atom[ito].pos[o2], cel.atom[ito].pos[o1]);
			// proton transfer coordinate
			double ptcoord = dis2 - dis3;

			// extremely important !!! mohan added 2017-01-18
			bool should_count=true;
			if(INPUT.system=="hydronium")
			{
				if(cel.snapshot_index < s1 and ptcoord > 0)
				{
					should_count=false;
				}
				if(cel.snapshot_index > s2 and ptcoord < 0)
				{
					should_count=false;
				}
			}
			else if(INPUT.system=="hydroxide")
			{
				if(cel.snapshot_index < s1 and ptcoord < 0)
				{
					should_count=false;
				}
				if(cel.snapshot_index > s2 and ptcoord > 0)
				{
					should_count=false;
				}
			}
			//				cout << setw(15) << dis1 << setw(15) << dis2 << setw(15) << dis3 << endl;

			int indexX;
			int indexY; 

			if(INPUT.func==1)
			{
				indexX = (ptcoord-INPUT.x0)/INPUT.dx;
				indexY = (dis1-INPUT.y0)/INPUT.dy;
			}
			else if(INPUT.func==0)
			{
				indexX = (ptcoord-INPUT.x0)/INPUT.dx;
				if(ptcoord<0.0)
					indexY = (angle1-INPUT.y0)/INPUT.dy;
				else if(ptcoord>0.0)
					indexY = (angle2-INPUT.y0)/INPUT.dy;
			}
			// calculate how the coordination number change during the proton transfer 
			else if(INPUT.func==-1)
			{
				assert(INPUT.ps_dt>0.0);
				double t1 = cel.snapshot_time;
				double t2 = PT.snapshot_time_pt[ip-1];
				int iii = (t1-t2)/INPUT.ps_dt+INPUT.npoints;
				if(iii>=0 and iii<this->np_tot)	
				{
					if(INPUT.system=="hydroxide")
					{
						if(water[o1].nH==1)
						{
							coor_num[iii] += water[o1].naccept;	
							pt_time[iii] += 1.0;
						}
						else if(water[o2].nH==1)
						{
							coor_num[iii] += water[o2].naccept;	
							pt_time[iii] += 1.0;
						}
					}
				}
			}

			//cout << setw(10) << indexX << setw(10) << indexY << endl;

			if(indexX<INPUT.nx and indexY<INPUT.ny and indexX>=0 and indexY>=0 and should_count==true)
			{
				coord_xy[indexX][indexY]+=1.0;
			}
		}
	}
	return;
}
