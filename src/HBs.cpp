#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "math.h"

HBs::HBs()
{
	nacc = -1;
	ndon = -1;
	last_nacc = -1;
	last_ndon = -1;

	snapshot_index = -1;
	snapshot_time = -1.0;

	last_snapshot_index = -1;
	last_snapshot_time = -1.0;

	last_nH=-1;
	last_Hindex=new int[8];
}

HBs::~HBs()
{
	delete[] last_Hindex;
}

void HBs::Routine()
{
	TITLE("HBs","Routine");
	
	ofs_running << "Compute the H-bond Network Information of Water" << endl;

	//************
	// OPEN FILES 
	//************
	ofs_bond.open("bond.dat");
	ofs_hbcase.open("hbcase.dat");

	ofs_en.open("en.dat");

	ofs_trans.open("trans.dat");
	ofs_wire.open("waterwire.dat");
	ofs_zundel.open("zundel.dat"); // compute how many zundel-like configurations
	ofstream ofs_angle("solvation_angle.dat");
	ofstream ofs_sd("solvation_distance.dat");
	ofs_running << "open file bond.dat" << endl;
	ofs_running << "open file hbcase.dat" << endl;
	ofs_running << "open file trans.dat" << endl;
	ofs_running << "open file waterwire.dat" << endl;
	ofs_running << "open file zundel.dat" << endl;
	ofs_running << "open file solvation_angle.dat" << endl;
	ofs_running << "open file solvation_distance.dat" << endl;

	this->avg_aion = 0;
	this->avg_dion = 0;
	this->avg_count_ion = 0; // counting for ions

	this->avg_awater = 0;
	this->avg_dwater = 0;
	this->avg_count_water = 0; // counting for water molecules

	this->accept_count = new int[10]();
	this->donate_count = new int[10]();
	this->HB_count = new int[10]();
	this->HB_count_acc = new int*[10];
	for (int ih=0; ih<10; ih++)
	{
		HB_count_acc[ih] = new int[10];
		for (int ih2 = 0; ih2<10; ih2++)
		{
			HB_count_acc[ih][ih2] = 0;
		}
	}

	this->old_ion_index=-1;

	// setup for angles	
	assert(INPUT.d_angle>0);
	const int nangle = INPUT.acut_hoo/INPUT.d_angle+10;
	assert(nangle>0);
	angle_ion = new double[nangle]();
	angle_water = new double[nangle];

	// setup for distances
	assert(INPUT.dr>0);
	const int nr = INPUT.rcut_oo/INPUT.dr+10;
	assert(nr>0);
	accept_dis = new double[nr];
	donate_dis = new double[nr];
	for(int ir=0; ir<nr; ++ir) accept_dis[ir]=0.0;
	for(int ir=0; ir<nr; ++ir) donate_dis[ir]=0.0;
	
	//**************************
	// BEGIN CALCULATING DATA
	//**************************
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

		if(cel.read_and_used==false) 
		{
			cel.clean(); // renxi added 20200614
			continue;
		}
		++count_geometry_number;
		cout << "snapshot " << igeo << endl;
		search_hydrogen_bonds(cel, igeo);
	 	cel.clean();	
	}	

	//********************
	//	PRINT INFORMATION
	//********************

	// clean up those who have finshed their jobs	
	ofs_bond.close();
	ofs_hbcase.close();
	ofs_en.close();
	ofs_trans.close();
	ofs_wire.close();
	ofs_zundel.close();

	//--------------------------------------------------------------------------
	// print out the angles of H-O-O for solvation structure of water and ion 
	//--------------------------------------------------------------------------
	ofs_angle << "Degrees Ion Water" << endl;
	double sum1=0.0;
	double sum2=0.0;
	for(int i=0; i<nangle; ++i)
	{
		sum1+=angle_ion[i];
		sum2+=angle_water[i];
	}
	sum1*=INPUT.d_angle;
	sum2*=INPUT.d_angle;
	//cout << sum1 << " " << sum2 << endl;
	if(sum1>0.0 and sum2>0.0)
	{
		for(int i=0; i<nangle; ++i)
		{
			double now_angle = i*INPUT.d_angle;
			angle_ion[i]/=sum1;
			angle_water[i]/=sum2;
			if(now_angle <= INPUT.acut_hoo )
			{		
				ofs_angle << now_angle << " " << angle_ion[i] << " " << angle_water[i] << endl;
			}
		}
	}

	// clean up
	//cout << "delete angle_ion and angle_water" << endl;
	//delete[] angle_ion;
	//delete[] angle_water;
	ofs_angle.close();
	//cout << "2" << endl;

	//------------------------------
	// accept_dis and donate_dis 
	//------------------------------
	ofs_sd << "index accept_dis donate_dis" << endl;
	double suma=0.0;
	double sumd=0.0;
	for(int i=0; i<nr; ++i)
	{	
		suma+=accept_dis[i];
		sumd+=donate_dis[i];
	}
	suma*=INPUT.dr;
	sumd*=INPUT.dr;
	if(suma>0 and sumd>0)
	{
		for(int i=0; i<nr; ++i)
		{
			double r = INPUT.dr * i;
			if( r <= INPUT.rcut_oo )
			{
				accept_dis[i]/=suma;
				donate_dis[i]/=sumd;
				ofs_sd << r << " " << accept_dis[i] << " " << donate_dis[i] << endl;
			}
		}
	}

	// clean up
	delete[] accept_dis;
	delete[] donate_dis;
	ofs_sd.close();

	//------------------------------
	// print out data for ions
	//------------------------------
	ofs_running << " --- THE H-BOND NETWORK OF IONS ---" << endl;
	suma=0.0;
	double sumb=0.0;
	for(int i=0; i<10; ++i) suma += accept_count[i];
	for(int i=0; i<10; ++i) sumb += donate_count[i];
	if(suma>0.0)
	{	
		ofs_running << setw(10) << "HB_number"  << setw(20) << "accept_count" 
					<< setw(20) << "percentage(%)" << endl;
		for(int i=0; i<10; ++i)
		{
			ofs_running << setw(10) << i << setw(20) << accept_count[i] 
					<< setw(20) << accept_count[i]/suma*100 << endl;
		}
	}
	if(sumb>0.0)
	{	
		ofs_running << setw(10) << "HB_number"  << setw(20) << "donate_count" 
					<< setw(20) << "percentage(%)" << endl;
		for(int i=0; i<10; ++i)
		{
			ofs_running << setw(10) << i << setw(20) << donate_count[i] 
					<< setw(20) << donate_count[i]/sumb*100 << endl;
		}
	}

	if(avg_count_ion>0)
	{
		ofs_running << " Counting for ions (nH=1 or 3) = " << avg_count_ion << endl;
		ofs_running << " Average accepted HBs of ions " << avg_aion/avg_count_ion << endl;
		ofs_running << " Average donating HBs of ions " << avg_dion/avg_count_ion << endl;
	}
	
	//----------------------------------
	// print out data for water systems
	//----------------------------------
	ofs_running << " --- THE H-BOND NETWORK OF WATER MOLECLUES ---" << endl;
	double sum=0.0;
	for(int i=0; i<10; ++i)
	{
		sum += HB_count[i];
	}

	ofstream ofs_HBacc("HB_acc.txt");
	if(sum>0.0)
	{	
		ofs_running << setw(10) << "HB_number" << setw(20) << "HB_count" << setw(20) << "percentage(%)" << endl;
		ofs_HBacc << setw(10) << "HB_number" << setw(20) << "HB_count" << setw(20) << "percentage(%)" ;
		for (int iacc=0; iacc<10; iacc++)
		{
			ofs_HBacc << setw(20) << "acc = " << iacc;
		}
		ofs_HBacc << endl;
		for(int i=0; i<10; ++i)
		{
			ofs_running << setw(10) << i << setw(20) << HB_count[i] << setw(20) << HB_count[i]/sum*100 << endl;
			ofs_HBacc << setw(10) << i << setw(20) << HB_count[i] << setw(20) << HB_count[i]/sum*100;
			for(int iacc=0; iacc<10; iacc++)
			{
				if (HB_count[i] > 0)
				{
					ofs_HBacc << setw(20) << (double)(HB_count_acc[i][iacc])/(double)(HB_count[i]) * 100 ;
				}
				else
				{
					ofs_HBacc << setw(20) << 0 ;
				}
			}
			ofs_HBacc << endl;
		}
		ofs_HBacc.close();

	}

	if(avg_count_water>0)
	{
		ofs_running << " Average accept HBs for water molecules is " << avg_awater/avg_count_water << endl;
		ofs_running << " Average donate HBs for water molecules is " << avg_dwater/avg_count_water << endl;
	}


	delete[] accept_count;
	delete[] donate_count;
	delete[] HB_count;
	for (int ih=0; ih<10; ih++)
	{
		delete[] HB_count_acc[ih];
	}
	delete[] HB_count_acc;

	return;
}


void HBs::search_hydrogen_bonds(const Cell &cel, const int &igeo)
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
	//----------------------------------------------------------------------------------------------------------------------
	// print out waterwires that ions are involved in
	//----------------------------------------------------------------------------------------------------------------------
	for(int o1=0; o1<cel.atom[ito].na; ++o1)
	{
		if(INPUT.system=="hydroxide" and water[o1].nH==1 and Water::nions==1)
		{
			ofs_wire << cel.snapshot_index << " " << cel.snapshot_time << endl;
			ofs_running << "ion_position " << cel.snapshot_time << " " << cel.atom[ito].pos[o1].z << endl; // mohan added
			int count=0;
			for(int iacc2=0; iacc2<water[o1].naccept; ++iacc2)
			{
				int o2 = water[o1].acceptO[iacc2];
				int h1 = water[o1].acceptH[iacc2];
				if(water[o2].nH!=2) continue;
				for(int iacc3=0; iacc3<water[o2].naccept; ++iacc3)
				{
					int o3 = water[o2].acceptO[iacc3];
					int h2 = water[o2].acceptH[iacc3];
					if(water[o3].nH!=2) continue;
					if(o3==o1) continue;
					++count;
				}
			}
			ofs_wire << count << endl;
			for(int iacc2=0; iacc2<water[o1].naccept; ++iacc2)
			{
				int o2 = water[o1].acceptO[iacc2];
				int h1 = water[o1].acceptH[iacc2];
				if(water[o2].nH!=2) continue;
				for(int iacc3=0; iacc3<water[o2].naccept; ++iacc3)
				{
					int o3 = water[o2].acceptO[iacc3];
					int h2 = water[o2].acceptH[iacc3];
					if(water[o3].nH!=2) continue;
					if(o3==o1) continue;
					ofs_wire << "O123 " << o1+1 << " " << o2+1 << " " << o3+1 << " h12 " << h1+1 << " " << h2+1 << endl;
				}
			}
		}
		else if(INPUT.system=="hydronium" and water[o1].nH==3 and Water::nions==1) 
		{	
			ofs_wire << cel.snapshot_index << " " << cel.snapshot_time << endl;
			int count=0;
			for(int iacc2=0; iacc2<water[o1].ndonate; ++iacc2)
			{
				int o2 = water[o1].donateO[iacc2];
				int h1 = water[o1].donateH[iacc2];
				if(water[o2].nH!=2) continue;
				for(int iacc3=0; iacc3<water[o2].ndonate; ++iacc3)
				{
					int o3 = water[o2].donateO[iacc3];
					int h2 = water[o2].donateH[iacc3];
					if(water[o3].nH!=2) continue;
					if(o3==o1) continue;
					++count;
				}
			}
			ofs_wire << count << endl;
			for(int iacc2=0; iacc2<water[o1].ndonate; ++iacc2)
			{
				int o2 = water[o1].donateO[iacc2];
				int h1 = water[o1].donateH[iacc2];
				if(water[o2].nH!=2) continue;
				for(int iacc3=0; iacc3<water[o2].ndonate; ++iacc3)
				{
					int o3 = water[o2].donateO[iacc3];
					int h2 = water[o2].donateH[iacc3];
					if(water[o3].nH!=2) continue;
					if(o3==o1) continue;
					ofs_wire << "O123 " << o1+1 << " " << o2+1 << " " << o3+1 << " h12 " << h1+1 << " " << h2+1 << endl;
				}
			}
		}
	}
	// compute zundel-like configurations
	int o2_save = -1;
	double dis_delta_save = 100.0;
	for(int o1=0; o1<cel.atom[ito].na; ++o1)
	{
		if(INPUT.system=="hydronium" and water[o1].nH==3)
		{ 
			for(int iacc2=0; iacc2<water[o1].ndonate; ++iacc2)
			{
				int o2 = water[o1].donateO[iacc2];
				int h1 = water[o1].donateH[iacc2];
			    double dis_o1_h1 = distance(cel.atom[ito].pos[o1], cel.atom[ith].pos[h1], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
			    double dis_o2_h1 = distance(cel.atom[ito].pos[o2], cel.atom[ith].pos[h1], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				double dis_delta = abs(dis_o1_h1 - dis_o2_h1);
			//	ofs_zundel << setw(12) << cel.snapshot_index << 
			//	setw(5) << o2 << setw(20) << dis_o1_h1 << setw(20) << dis_o2_h1 << setw(20) << dis_delta << endl; 

				if(dis_delta<dis_delta_save)
				{
					o2_save = o2;
					dis_delta_save = dis_delta;
				}
			}
		}
	}	
	ofs_zundel << setw(12) << cel.snapshot_index << 
		setw(5) << o2_save << setw(20) << dis_delta_save << endl; 

	//---------------------------------
	// ofs_bond, ofs_hbcase, ofs_trans
	//---------------------------------
	ofs_bond << cel.snapshot_index << " " << cel.snapshot_time << endl;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		ofs_bond << " atom" << setw(5) << ia+1 << " A" << setw(2) << water[ia].naccept 
		<< " D" << setw(2) << water[ia].ndonate;

		for(int i=0; i<8; ++i) ofs_bond << setw(5) << water[ia].acceptO[i]+1;
		for(int i=0; i<4; ++i) ofs_bond << setw(5) << water[ia].donateO[i]+1;
		ofs_bond << setw(10) << cel.snapshot_index << " " << cel.snapshot_time << endl;

		if( (INPUT.system=="hydroxide" and water[ia].nH==1 and Water::nions==1) or
            (INPUT.system=="hydronium" and water[ia].nH==3 and Water::nions==1) )
		{
			// nacc of last snapshot of previous ion
			last_nacc = nacc; // previous number of accepted HBs for ion
			last_ndon = ndon;
			last_snapshot_index = snapshot_index;
			last_snapshot_time = snapshot_time; 

			nacc = water[ia].naccept; // ia is now ion
			ndon = water[ia].ndonate; // 
			snapshot_index = cel.snapshot_index;
			snapshot_time = cel.snapshot_time;
	
			// mohan add 2017-1-2
			// need to know which proton is involved in the proton transfer event
			int h_index=-1;
			if(INPUT.system=="hydronium")
			{
				bool foundH = false;
				for(int ih=0; ih<water[ia].nH; ++ih)
				{
					if(foundH) break;
					for(int ih2=0; ih2<last_nH; ++ih2)
					{
						//cout << setw(10) << water[ia].indexH[ih] << setw(10) << last_Hindex[ih2] << endl;
						if(water[ia].indexH[ih]==last_Hindex[ih2])
						{
							h_index = water[ia].indexH[ih];
							foundH=true;
							break;
						}
					}
				}
			} // end hydronium
			else if(INPUT.system=="hydroxide")
			{
				if(old_ion_index!=-1 and ia!=old_ion_index)
				{
					if(water[old_ion_index].nH!=2)
					{
						cout << "snapshot " << cel.snapshot_index << endl;
						cout << "ia=" << ia << endl;
						cout << "old_ion_index=" << old_ion_index << endl;
						cout << "water[old_ion_index].nH=" << water[old_ion_index].nH << endl;
					}
					for(int ih=0; ih<water[old_ion_index].nH;++ih)
					{
						for(int ih2=0; ih2<last_nH; ++ih2)
						{
							if(water[old_ion_index].indexH[ih]!=last_Hindex[0])
							{
								h_index = water[old_ion_index].indexH[ih];
								break;
							}
						}
					}
				}
			} // end hydroxide

			assert(water[ia].nH<=8);
			last_nH = water[ia].nH;
			for(int ih=0; ih<water[ia].nH; ++ih)
			{
				last_Hindex[ih]=water[ia].indexH[ih];
			}

			if(ia!=old_ion_index)
			{
				if(old_ion_index!=-1)
				{
					ofs_trans << setprecision(12) << setw(10) << this->snapshot_index 
						<< setw(15) << this->snapshot_time
						<< setw(5) << ia+1 // ion index starting from 1, not 0 
						<< setw(3) << nacc 
						<< setw(3) << ndon 
						<< setw(10) << last_snapshot_index
						<< setw(15) << last_snapshot_time
						<< setw(5) << old_ion_index+1
						<< setw(3) << last_nacc 
						<< setw(3) << last_ndon 	
						<< setw(5) << h_index+1
						<< endl;
					//assert(h_index>=0);
				}
				old_ion_index=ia;
			}

			assert(nacc<10);
			assert(ndon<10);
			++accept_count[nacc]; 
			++donate_count[ndon];

			ofs_hbcase << setprecision(12) << setw(10) << count_geometry_number << setw(10) << cel.snapshot_index << setw(15) << cel.snapshot_time
				<< setw(5) << ia+1 << setw(3) << nacc  << setw(3) << ndon;

			for(int iii=0; iii<8; ++iii) ofs_hbcase << setw(5) << water[ia].acceptO[iii]+1;
			for(int iii=0; iii<4; ++iii) ofs_hbcase << setw(5) << water[ia].donateO[iii]+1;
			ofs_hbcase << endl;
		} // for water
		else
		{
			const int total_HB = water[ia].naccept+water[ia].ndonate;
			++HB_count[total_HB];
			++HB_count_acc[total_HB][water[ia].naccept];
		}
	}


	// output information
	if(INPUT.system == "hydroxide" or INPUT.system == "hydronium")
	{
		for(int ia=0; ia<cel.atom[ito].na; ++ia)
		{
			// if the system has hydronium or hydroxide
			if( (water[ia].nH==1 and INPUT.system=="hydroxide" and Water::nions==1) or 
                (water[ia].nH==3 and INPUT.system=="hydronium" and Water::nions==1) )
// renxi changed 20210110
//			if( (water[ia].nH==1 and INPUT.system=="hydroxide" ) or 
//                (water[ia].nH==3 and INPUT.system=="hydronium" ) )
			{
				// the default value of INPUT.nacc is -1.
				// if the INPUT.nacc !=1, only print out information for ions who have
				// 'INPUT.nacc' accetped H bonds.
				if(INPUT.nacc==-1 or (INPUT.nacc!=-1 and INPUT.nacc==water[ia].naccept) )
				{
					for(int iacc=0; iacc < water[ia].naccept; ++iacc)
					{
						double angle0=water[ia].accept_angle[iacc];
						int iangle = (int)(angle0/INPUT.d_angle);
						angle_ion[iangle]++;
						double dis0 = water[ia].accept_disO[iacc];
						int iaccept = (int)(dis0/INPUT.dr);
						if(dis0<INPUT.rcut_oo) accept_dis[iaccept]++;
					}
					for(int idon=0; idon < water[ia].ndonate; ++idon)
					{
						double angle0=water[ia].donate_angle[idon];
						int iangle = (int)(angle0/INPUT.d_angle);
						angle_ion[iangle]++;
						double dis0 = water[ia].donate_disO[idon];
						int idonate = (int)(dis0/INPUT.dr);
						if(dis0<INPUT.rcut_oo) donate_dis[idonate]++;
					}

					avg_aion+=water[ia].naccept;
					avg_dion+=water[ia].ndonate;
					avg_count_ion++;
				} // end nacc
			} // end nH
		} // end ia
	} // end hydronium and hydroxide

	// for all oxygen atoms
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		if(water[ia].nH==2)
		{
			for(int iacc=0; iacc < water[ia].naccept; ++iacc)
			{
				double angle0=water[ia].accept_angle[iacc];
				if(angle0 <= INPUT.acut_hoo)
				{
					int iangle = (int)(angle0/INPUT.d_angle);
					angle_water[iangle]++;
				}
			}
			for(int idon=0; idon < water[ia].ndonate; ++idon)
			{
				double angle0=water[ia].donate_angle[idon];
				if(angle0 <= INPUT.acut_hoo)
				{
					int iangle = (int)(angle0/INPUT.d_angle);
					angle_water[iangle]++;
				}
			}
			avg_awater+=water[ia].naccept;
			avg_dwater+=water[ia].ndonate;
			avg_count_water++;
		}
	}

	// record en information
	// mohan added 2018-04-27
	ofs_en << cel.snapshot_index << " " << cel.snapshot_time << endl;
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		if(water[ia].nH==2)
		{
			ofs_en << " O" << setw(5) << ia+1 << " H ";
			double sum_dis=0.0;
			for(int ih=0; ih<water[ia].nH; ++ih)
			{
				const int ind_h = water[ia].indexH[ih];
				double dis = distance(cel.atom[ito].pos[ia], cel.atom[ith].pos[ind_h], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
				ofs_en << setw(5) << ind_h+1;
				ofs_en << setw(15) << dis;
				sum_dis+=dis;
			}
			ofs_en << setw(15) << sum_dis;
			ofs_en << endl;
		}
	}
	delete[] water;
	return;
}

bool HBs::accepted(const Cell &cel, Water *water, const int &ito, const int &ith, 
	const int &ia, const int &ia2, double &angle0, int &H_id, double &dis)
{
	for(int ih=0; ih<water[ia2].nH; ++ih)
	{
		const int ind = water[ia2].indexH[ih];
		// angle formed by H-O2-O1
		angle0 = angle(cel, cel.atom[ith].pos[ind], cel.atom[ito].pos[ia2], cel.atom[ito].pos[ia]);
		if (INPUT.HB_defination==1)
		{
			if( angle0 < INPUT.acut_hoo ) // unit is degrees
			{
				H_id = ind;
				return true;
			}
		}
		else if (INPUT.HB_defination==2)
		{
			assert(INPUT.HB_defination_k != 1000 and INPUT.HB_defination_b != -1);
			if (angle0 < INPUT.HB_defination_k*dis+INPUT.HB_defination_b)
			{
				H_id = ind;
				return true;
			}
		}
	}
	return false;
}

bool HBs::donate(const Cell &cel, Water *water, const int &ito, const int &ith, const int &ia, 
	const int &ia2, double &angle0, int &H_id, double &dis)
{
	for(int ih=0; ih<water[ia].nH; ++ih)
	{
		const int ind = water[ia].indexH[ih];
		// angle formed by H-O2-O1
		angle0 = angle(cel, cel.atom[ith].pos[ind], cel.atom[ito].pos[ia], cel.atom[ito].pos[ia2]);
		if (INPUT.HB_defination==1)
        {
			if( angle0 < INPUT.acut_hoo ) // unit is degrees
			{
				H_id = ind;
				return true;
			}
		}
		else if (INPUT.HB_defination == 2)
		{
                        assert(INPUT.HB_defination_k != 1000 and INPUT.HB_defination_b != -1);
                        if (angle0 < INPUT.HB_defination_k*dis+INPUT.HB_defination_b)
                        {
                                H_id = ind;
                                return true;
                        }
                }
	}
	return false;
}


// center atom is pos2
double HBs::angle(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &pos3) 
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


void HBs::setup_water(const Cell &cel, Water *water)
{
	// obtain ito and ith
	int ito=-1;
	int ith=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
	}
	assert(ito>=0); assert(ith>=0);

	//cout << "celldm123" << endl;
	//cout << INPUT.celldm1 << endl;
	//cout << INPUT.celldm2 << endl;
	//cout << INPUT.celldm3 << endl;

	// search for O-H pair

	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		//cout << cel.atom[ito].pos[ia].x << " " <<  cel.atom[ito].pos[ia].y << " " << cel.atom[ito].pos[ia].z << endl;
		for(int ia2=0; ia2<cel.atom[ith].na; ++ia2)
		{
			double dis = distance(cel.atom[ito].pos[ia], cel.atom[ith].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
			//cout << "distance of " << ia2 << " is " << dis << endl;
			if(dis < INPUT.rcut_oh)
			{
				int ind=water[ia].nH;
				water[ia].indexH[ind]=ia2;
				water[ia].disH[ind]=dis;
				water[ia].nH++;	
			}
		}
		//	cout << water[ia].nH << endl;
		//cout << "water: " << ia << " number of H: " << water[ia].nH << endl;

		if( (INPUT.system=="hydroxide" and water[ia].nH==1) or
            (INPUT.system=="hydronium" and water[ia].nH==3) ) 
		{
			Water::nions++;

			// find out how many O around H
			if(INPUT.system=="hydronium")
			{
				for(int ih=0; ih<water[ia].nH; ++ih)
				{
					const int iah = water[ia].indexH[ih];
					for(int ia3=0; ia3<cel.atom[ito].na; ++ia3)
					{
						double dis = distance(cel.atom[ito].pos[ia3], cel.atom[ith].pos[iah], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						if(dis<1.35) water[ia].nO_H++;
					}
				}
//				ofs_running << "number_of_O is " << water[ia].nO_H << endl;
			}
		}
		//cout << ia << " " << water[ia].indexH[0] << " " << water[ia].indexH[1] << endl;

	}


	//cout << "First water has " << water[0].nH << " H atoms." << endl;
	// search for O-O pair
	for(int ia=0; ia<cel.atom[ito].na; ++ia)
	{
		for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
		{
			if(ia==ia2) continue;
			double dis = distance(cel.atom[ito].pos[ia], cel.atom[ito].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
			if(dis < INPUT.rcut_oo)
			{
				water[ia].nO++;
				double angle0=-1.0;
				int tmp_idH=-1;
				if( accepted(cel, water, ito, ith, ia, ia2, angle0, tmp_idH, dis) )
				{	
					int ind = water[ia].naccept;
					water[ia].acceptO[ind] = ia2;
					water[ia].acceptH[ind] = tmp_idH;
					water[ia].accept_angle[ind] = angle0;
					water[ia].accept_disO[ind] = dis;
					water[ia].naccept++;
				}
				tmp_idH=-1;
				if( donate(cel, water, ito, ith, ia, ia2, angle0, tmp_idH, dis) )
				{	
					int ind = water[ia].ndonate;
					water[ia].donateO[ind] = ia2;
					water[ia].donateH[ind] = tmp_idH;
					water[ia].donate_angle[ind] = angle0;
					water[ia].donate_disO[ind] = dis;
					water[ia].ndonate++;
				}
				//cout << water[ia2].indexH[0] << " " << water[ia2].indexH[1] << endl;
				tmp_idH = -1;
			}
		}
	}

	return;
}
