#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "HBs_near_PT.h"
#include <math.h>

HBs_near_PT::HBs_near_PT(){

}

HBs_near_PT::~HBs_near_PT(){

}

void HBs_near_PT::Routine(){
	//assert(INPUT.system=="hydroxide");
	//ifstream ifs_trans("trans.dat");
	//ifstream ifs_bond("bond.dat");
	cout << INPUT.func_b << endl;
	ifstream ifs_trans("trans.dat");
	ifstream ifs_mj("mj.dat");
	ofstream ofs("HBs_near_PT.txt");
	int n = INPUT.msd_dt0/INPUT.msd_dt; 
	string useless;
	int is_doc = 0;
	hbs_bf_pt = new double[n+1](); // store HBs before PT;
	hbs_aft_pt = new double[n+1]();
	count_bf_pt = new int[n+1]();
	count_aft_pt = new int[n+1]();
	//hbs_total = new double[2*n+2]();

	snapshot_pt = new int[INPUT.nbands]();
	snapshot_time_pt = new double[INPUT.nbands]();
	iindex_p = new int[INPUT.nbands]();
	iindex = new int[INPUT.nbands]();
	return_jump = new bool[INPUT.nbands];
	int orig_iidx;
	ifs_trans >> useless >> useless >> useless >> useless >> useless >> useless >> useless >> orig_iidx >> useless >> useless >> useless;
	orig_iidx--;
	for(int is = 0; is<INPUT.nbands; is++)
	{
		int ii, ii_p;
		int ss_pt;
		double ss_time_pt;
		string pt_type;
		ifs_trans >> useless >> useless >> ii >> useless >> useless >> ss_pt >> ss_time_pt >> ii_p >> useless >> useless >> useless;
		for(int i = 0; i<12; i++){ifs_mj >> useless;}
		ifs_mj >> pt_type;
		cout << pt_type << endl;
		//if(pt_type=="rattling"){continue;}
		if(INPUT.func==2 and pt_type != "single")
		{
			continue;
		}
		else if(INPUT.func==3 and pt_type != "double")
		{
			continue;
		}
		else if(INPUT.func==4 and pt_type != "triple")
		{
			continue;
		}
		snapshot_pt[is_doc] = ss_pt;
		snapshot_time_pt[is_doc] = ss_time_pt;
		iindex[is_doc] = ii-1;
		iindex_p[is_doc] = ii_p-1;
		//cout << is_doc << endl;
		cout << iindex[is_doc] << " " << iindex_p[is_doc] << " " << snapshot_pt[is_doc] << " " << snapshot_time_pt[is_doc] << endl;
		is_doc++;
		
	}

	int nreturn_jump = 0;
	//cout << "iindex[0] = " << iindex[0] << endl;
	//cout << "orig_iidx = " << orig_iidx << endl;
	if (iindex[0] == orig_iidx) 
	{
		this->return_jump[0] = true;
		nreturn_jump++;
		cout << 0 << " " << true << endl;
	}
	else 
	{
		this->return_jump[0] = false;
		cout << 0 << " " << false << endl;
	}
	for (int is=1; is<INPUT.nbands; is++)
	{
		if (iindex[is] == iindex_p[is-1])
		{
		 	this->return_jump[is] = true;
			nreturn_jump++;
			cout << is << " " << true << endl;
		}
		else 
		{
			this->return_jump[is] = false;
			cout << is << " " << false << endl;
		}
	}
	INPUT.nbands = is_doc;
	/*
	if (INPUT.non_return > 0)
	{
		INPUT.nbands -= nreturn_jump;
	}
	*/
	cout << INPUT.nbands << endl;

	for(int i = 0; i < n+1; i++)
	{
		hbs_aft_pt[i] = 0;
		hbs_bf_pt[i] = 0;
		count_bf_pt[i] = 0;
		count_aft_pt[i] = 0;
	}
	//for(int i = 0; i < 2*n+2; i++)
	//{
	//	hbs_total[i] = 0;
	//}
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
		for(int is = 0; is<INPUT.nbands; is++)
		{
			if (INPUT.non_return > 0 and this->return_jump[is] == true)
			{
				continue;
			}
			if(cel.snapshot_time <= snapshot_time_pt[is] and cel.snapshot_time > snapshot_time_pt[is]-n*INPUT.msd_dt)
			{
				int ig = round((snapshot_time_pt[is]-cel.snapshot_time)/INPUT.msd_dt);
				if((water[iindex_p[is]].nH==1 and INPUT.system == "hydroxide") or (water[iindex_p[is]].nH==3 and INPUT.system == "hydronium"))
				{
					if(INPUT.func_b == 2)// record new ion
					{
						hbs_bf_pt[n-ig] += water[iindex[is]].naccept;
						count_bf_pt[n-ig]++;
					}
					else
					{
						hbs_bf_pt[n-ig] += water[iindex_p[is]].naccept;
						count_bf_pt[n-ig]++;
					}
				//cout << "1 " << snapshot_time_pt[is]-n*INPUT.msd_dt << " " << cel.snapshot_time << " " << snapshot_time_pt[is] << " " << water[iindex_p[is]].naccept << " " << water[iindex[is]].naccept << endl;
				}
				//hbs_bf_pt[n-ig] += water[iindex_p[is]].naccept;
				//count_bf_pt[n-ig]++;
			}
			//cout << 1 << endl;
			else if(cel.snapshot_time > snapshot_time_pt[is] and cel.snapshot_time < snapshot_time_pt[is]+n*INPUT.msd_dt)
			{
				int ig = round((cel.snapshot_time-snapshot_time_pt[is])/INPUT.msd_dt);
				if((water[iindex[is]].nH==1 and INPUT.system == "hydroxide") or (water[iindex[is]].nH==3 and INPUT.system == "hydronium" ))
				{
					if(INPUT.func_b == 3)// record old ion
					{
						hbs_aft_pt[ig-1] += water[iindex_p[is]].naccept;
						count_aft_pt[ig-1]++;
					}
					else
					{
						hbs_aft_pt[ig-1] += water[iindex[is]].naccept;
						count_aft_pt[ig-1]++;
					}
					
				}
				//hbs_aft_pt[ig-1] += water[iindex[is]].naccept;
				//count_aft_pt[ig-1]++;
				//cout << "2 " << snapshot_time_pt[is] << " " << cel.snapshot_time << " " << snapshot_time_pt[is]+n*INPUT.msd_dt << " " << water[iindex_p[is]].naccept << " " << water[iindex[is]].naccept << endl;
			}

		}
		delete[] water;
	}// for igeo
//	for(int i=0; i<INPUT.nbands; i++)
//	{
//		cout << iindex[i] << " " << iindex_p[i] << endl;
//	}
	for(int i=0; i<=n; i++)
	{
		cout << INPUT.msd_dt*(i-n) << " " << count_bf_pt[i] << " " << hbs_bf_pt[i] << endl;
	}
	for(int i=0; i<n; i++)
	{
		cout << INPUT.msd_dt*(i+1) << " " << count_aft_pt[i] << " " << hbs_aft_pt[i] << endl;
	}

	for(int i=0; i<=n; i++)
	{
		if(count_bf_pt[i]<0 or hbs_bf_pt[i]<0)
		{
			cout << "Warning! HBs before PT is wrong." << endl;
			cout << "i = " << i << endl;
			cout << "count_bf_pt[i] = " << count_bf_pt[i] << endl;
			cout << "hbs_bf_pt[i] = " << hbs_bf_pt[i] << endl;
			exit(0);
		}
		hbs_bf_pt[i] = hbs_bf_pt[i]/count_bf_pt[i];
		ofs << INPUT.msd_dt*(i-n) << " " << hbs_bf_pt[i] << endl;
	}
	for(int i=0; i<n; i++)
	{
		if(count_aft_pt[i]<0 or hbs_aft_pt[i]<0)
		{
			cout << "Warning! HBs after PT is wrong." << endl;
			cout << "i = " << i << endl;
			cout << "count_aft_pt[i] = " << count_aft_pt[i] << endl;
			cout << "hbs_aft_pt[i] = " << hbs_aft_pt[i] << endl;
			exit(0);
		}
		hbs_aft_pt[i] = hbs_aft_pt[i]/count_aft_pt[i];
		ofs << INPUT.msd_dt*(i+1) << " " << hbs_aft_pt[i] << endl;
	}
	ofs.close();
/*
	delete[] hbs_bf_pt;
	delete[] hbs_aft_pt;
	delete[] count_bf_pt;
	delete[] count_aft_pt;
	delete[] snapshot_time_pt;
	delete[] snapshot_pt;
	delete[] iindex;
	delete[] iindex_p;
	*/
}

int HBs_near_PT::round(const double &r)
{
    if(r - floor(r)>=0.5)
    {
    	return ceil(r);
    }
    else
    {
    	return floor(r);
    }
}
