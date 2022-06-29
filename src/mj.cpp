#include "cellFile.h"
#include "input.h"
#include "mj.h"
#include "math.h"

ProtonTransfer::ProtonTransfer()
{
	this->use_pt=false; 
	this->npt = 1;

	// proton transfer 
	snapshot_index = new int[npt];
	snapshot_time = new double[npt];	
	ion_index = new int[npt];
	nacc = new int[npt];
	ndon = new int[npt];

	// after a resting period; right before next proton transfer
	snapshot_index_pt = new int[npt];
	snapshot_time_pt = new double[npt];
	ion_index_pt = new int[npt];
	nacc_pt = new int[npt];
	ndon_pt = new int[npt];

	indexH = new int[npt];

	type_pt = new string[npt];
	ds_pt = new int[npt];
	dt_pt = new double[npt];
	ds_min = new int[npt];
	ds_max = new int[npt];
}

ProtonTransfer::~ProtonTransfer()
{
	delete[] snapshot_index;
	delete[] snapshot_time;
	delete[] ion_index;
	delete[] nacc;
	delete[] ndon;
	
	delete[] snapshot_index_pt;
	delete[] snapshot_time_pt;
	delete[] nacc_pt;
	delete[] ndon_pt;

	delete[] indexH;

	delete[] type_pt;
	delete[] ds_pt;
	delete[] dt_pt;
	delete[] ds_min;
	delete[] ds_max;
}

void ProtonTransfer::setup_PT()
{
	ifstream ifs_mj("mj.dat");
	if(!ifs_mj)
	{
		cout << "Cannot find the file: mj.dat, use_pt=false" << endl;
		use_pt=false;
		return;
	}
	else
	{
		use_pt=true;
	}

	this->npt=0;
    while(ifs_mj.eof()==0)
    {
		string tmp;
        READ_VALUE(ifs_mj, tmp);
        ++npt;
    }
	--npt;

	
	ofs_running << " Number of PT events recorded: " << npt << endl;
	ifs_mj.clear();
	ifs_mj.seekg(0);

	// initialize many arrays
	assert(npt>0);

	delete[] snapshot_index;
	delete[] snapshot_time;
	delete[] ion_index;
	delete[] nacc;
	delete[] ndon;

	delete[] snapshot_index_pt;
	delete[] snapshot_time_pt;
	delete[] ion_index_pt;
	delete[] nacc_pt;
	delete[] ndon_pt;

	delete[] indexH;

	delete[] type_pt;
	delete[] ds_pt;
	delete[] dt_pt;
	delete[] ds_min;
	delete[] ds_max;
	
	snapshot_index = new int[npt];
	snapshot_time = new double[npt];	
	ion_index = new int[npt];
	nacc = new int[npt];
	ndon = new int[npt];
	
	snapshot_index_pt = new int[npt];
	snapshot_time_pt = new double[npt];
	ion_index_pt = new int[npt];
	nacc_pt = new int[npt];
	ndon_pt = new int[npt];

	indexH = new int[npt];

	type_pt = new string[npt];
	ds_pt = new int[npt];
	dt_pt = new double[npt];
	ds_min = new int[npt];
	ds_max = new int[npt];

	// read in PT events
	int index=0;
	for(int i=0; i<npt; ++i)
	{
		ifs_mj >> index 
		>> snapshot_index[i]
		>> snapshot_time[i]
		>> ion_index[i]
		>> nacc[i]
		>> ndon[i]
		>> snapshot_index_pt[i]
		>> snapshot_time_pt[i]
		>> ion_index_pt[i]
		>> nacc_pt[i]
		>> ndon_pt[i]
		>> indexH[i]
		>> type_pt[i];
	//	cout << index << " " << type_pt[i] << endl;
	}


	// how much time it takes to transfer a proton from ion A to B
	for(int i=0; i<npt-1; ++i)
	{
		// how much time it takes to let proton transfer process occurs.
		dt_pt[i]=snapshot_time[i+1]-snapshot_time_pt[i];
		// how many snapshots it takes to let proton transfer occurs.
		ds_pt[i]=snapshot_index[i+1]-snapshot_index_pt[i];
		ds_min[i] = snapshot_index_pt[i]-ds_pt[i]*INPUT.ext_1-INPUT.movement_x;
		if(ds_min[i]<0) ds_min[i]=0;
		ds_max[i] = snapshot_index[i+1]+ds_pt[i]*INPUT.ext_2+INPUT.movement_y;
		ofs_running << "PT " << setw(5) << i+1 
		<< setw(10) << ds_pt[i] 
		<< setw(15) << dt_pt[i] 
		<< " dsmin " << setw(15) << ds_min[i]
		<< " dsmax " << setw(15) << ds_max[i]
		<< endl;
	}


	ifs_mj.close();
	return;
}


// find which PT event this "snapshot_i" belongs to
string ProtonTransfer::which_type_pt(const int &snapshot_i)
{
	string name_pt="none";
	for(int i=0; i<npt; ++i)
	{
		if(snapshot_i>=snapshot_index[i] and i==npt-1)
		{
			name_pt = type_pt[i];
		}
		else if(snapshot_i>=snapshot_index[i] and
			snapshot_i<snapshot_index[i+1])
		{
			name_pt = type_pt[i];
		}
		else
		{
			continue;
		}

		if(name_pt=="follow")
		{
			cout << "follow" << endl;
			int j=i;
			while(j>=0)
			{
				j--;
				if(type_pt[j]!="follow")
				{
					name_pt=type_pt[j];
					cout << "find " << name_pt << endl;
					break;
				}	
			}
		}

		if(name_pt!="none")
		{
			assert(name_pt!="follow");
			return name_pt;
		}
	}
}

MJ::MJ()
{
	dt_max=INPUT.dt_max;//unit is ps
}

MJ::~MJ()
{
}

void MJ::Routine()
{
	TITLE("MJ","Routine");
	
	cout << "Compute the multiple jumps information." << endl;

	ifstream ifs("trans.dat");

	this->ndim = INPUT.geo_2 - INPUT.geo_1;
	cout << "number of PT events: " << ndim << endl;	

	this->npt = new int[7]();
	this->npt_hyper = new int[7]();
	this->npt_tetra = new int[7]();

	int* si = new int[ndim]; // snapshot index
	double* st = new double[ndim]; // snapshot time
	int* ai = new int[ndim]; // ion index
	int* nacc = new int[ndim]; // number of accepting hydrogen bonds
	int* ndon = new int[ndim]; // number of donating hydrogen bonds

	int* last_si = new int[ndim]; // last snapshot index
	double* last_st = new double[ndim]; // last snapshot time
	int* last_ai = new int[ndim]; // last ion index
	int* last_nacc = new int[ndim]; // number of accepting hydrogen bonds
	int* last_ndon = new int[ndim]; // number of donating hydrogen bonds

	int* indexH = new int[ndim]; // index of H

	bool* done = new bool[ndim];
	
	for(int i=0; i<ndim; ++i)
	{
		ifs >> si[i] >> st[i] >> ai[i] >> nacc[i] >> ndon[i] 
		>> last_si[i] >> last_st[i] >> last_ai[i] >> last_nacc[i] >> last_ndon[i] >> indexH[i];

		done[i]=false;
	}

	ofstream ofs_mj("mj.dat");
	
	double total_time=0.0;
	for(int i1=0; i1<ndim-1; ++i1)
	{
		const int i2 = i1+1; // single jumps
		const int i3 = i2+1; // double jumps
		const int i4 = i3+1; // triple jumps
		const int i5 = i4+1; // quadraple jumps
		const int i6 = i5+1; // for quadraple jumps
		const int i7 = i6+1; // for five jumps
		const int i8 = i7+1; // for six jumps
		const int i9 = i8+1; // 

		bool sixj=false;
		bool fivej=false;
		bool quadraplej=false;
		bool triplej=false;
		bool doublej=false;
		bool singlej=false;
		bool rattling=false;

		ofs_mj << setprecision(12) << setw(5) << i1+1  // index for PT events
			<< setw(10) << si[i1]  // snapshot index
			<< setw(15) << st[i1]  // snapshot time
			<< setw(10) << ai[i1]  // index of ion (at time st[i1]) 
			<< setw(5) << nacc[i1] // number of accepted HB of this ion at time st[i1] 
			<< setw(5) << ndon[i1] // number of donating HB of this ion at time st[i1]
			<< setw(10) << last_si[i1+1]
			<< setw(15) << last_st[i1+1]
			<< setw(10) << last_ai[i1+1]
			<< setw(5) << last_nacc[i1+1] // the last snapshot of number of accepted HB before it leaves a1[i1]
			<< setw(5) << last_ndon[i1+1] // the last snapshot of number of donating HB before it leaves a1[i1] 
			<< setw(10) << indexH[i1] // index of proton involved in the proton transfer process
			<< setw(10); 
			// here one should use last_nacc from "i1+1"

		if(done[i1]==true) {ofs_mj << "follow" << endl; continue;}

		// six jump
		if(INPUT.func==1)
		{
			if(i7<ndim)
			{
				if(ai[i1] != ai[i3] and ai[i1] != ai[i4] and ai[i1] != ai[i5] and ai[i1] != ai[i6] and ai[i1] != ai[i7])
				{
					if(ai[i2] != ai[i4] and ai[i2] != ai[i5] and ai[i2] != ai[i6] and ai[i2] != ai[i7])
					{
						if(ai[i3] != ai[i5] and ai[i3] != ai[i6] and ai[i3] != ai[i7])
						{
							if(ai[i4] != ai[i6] and ai[i4] != ai[i7])
							{
								if(ai[i5] != ai[i7])
								{
									total_time = st[i7]-st[i2];
									if( total_time < dt_max )
									{
										cout << "possible six: " << si[i1] << " total time " << total_time << endl;;
										sixj=true;
										ofs_mj << "six" << endl;
										done[i1]=true;done[i2]=true;done[i3]=true;done[i4]=true;done[i5]=true;done[i6]=true;
										npt[6]++;
										if(INPUT.system == "hydroxide")
										{
											if(last_nacc[i1+1]>=4) npt_hyper[6]++; // pay special attention to this
											else if(last_nacc[i1+1]==3) npt_tetra[6]++;
										}
										else if(INPUT.system == "hydronium")
										{
											if(last_nacc[i1]==0) npt_hyper[6]++;
											else if(last_nacc[i1]==1) npt_tetra[6]++;
										}
									}
								}
							}
						}
					}
				}
			}
		}	
		if(sixj==true)continue;
		

		// five jumps
		if(INPUT.func==1)
		{
			if(sixj==false and i6<ndim)
			{
				if(ai[i1] != ai[i3] and ai[i1] != ai[i4] and ai[i1] != ai[i5] and ai[i1] != ai[i6])
				{
					if(ai[i2] != ai[i4] and ai[i2] != ai[i5] and ai[i2] != ai[i6])
					{
						if(ai[i3] != ai[i5] and ai[i3] != ai[i6])
						{
							if(ai[i4] != ai[i6])
							{
								total_time = st[i6]-st[i2];
								if( total_time < dt_max )
								{
									cout << "possible five: " << si[i1] << " total time " << total_time << endl;;
									fivej=true;
									ofs_mj << "five" << endl;
									done[i1]=true;done[i2]=true;done[i3]=true;done[i4]=true;done[i5]=true;
									npt[5]++;
									if(INPUT.system == "hydroxide")
									{
										if(last_nacc[i1+1]>=4) npt_hyper[5]++; // pay special attention to this
										else if(last_nacc[i1+1]==3) npt_tetra[5]++;
									}
									else if(INPUT.system == "hydronium")
									{
										if(last_nacc[i1]==0) npt_hyper[5]++; // pay special attention to this
										else if(last_nacc[i1]==1) npt_tetra[5]++;
									}
								}
							}
						}
					}
				}
			}
		}
		if(fivej==true)continue;


		// quadraple jump
		if(INPUT.func==1)
		{
			if(sixj==false and fivej==false and i5<ndim)
			{
				if(ai[i1] != ai[i3] and ai[i1] != ai[i4] and ai[i1] != ai[i5])
				{
					if(ai[i2] != ai[i4] and ai[i2] != ai[i5])
					{
						if(ai[i3] != ai[i5])
						{
							total_time = st[i5]-st[i2];
							if( total_time < dt_max )
							{
								cout << "possible quadraple: " << si[i1] << " total time " << total_time << endl;;
								quadraplej=true;
								ofs_mj << "quadraple" << endl;
								done[i1]=true;done[i2]=true;done[i3]=true;done[i4]=true;
								npt[4]++;
								if(INPUT.system == "hydroxide")
								{
									if(last_nacc[i1+1]>=4) npt_hyper[4]++; // pay special attention to this
									else if(last_nacc[i1+1]==3) npt_tetra[4]++;
								}
								else if(INPUT.system == "hydronium")
								{
									if(last_nacc[i1] == 0) npt_hyper[4]++;
									else if(last_nacc[i1] == 1) npt_tetra[4]++;
								}
							}
						}
					}
				}
			}
		}
		else if(INPUT.func==2)
		{
			if(sixj==false and fivej==false and i9<ndim)
			{
				if(ai[i1] != ai[i3] and ai[i1] != ai[i4] and ai[i1] != ai[i5])
				{
					if(ai[i2] != ai[i4] and ai[i2] != ai[i5])
					{
						if(ai[i3] != ai[i5])
						{
							total_time = st[i5]-st[i2];
							if( total_time < dt_max )
							{
								if(ai[i4]==ai[i6] and ai[i3]==ai[i7] and ai[i2]==ai[i8] and ai[i1]==ai[i9])
								{
									if( st[i9]-st[i5] > dt_max )
									{
										cout << "possible quadraple: " << si[i1] << " total time " << total_time << endl;;
										quadraplej=true;
										ofs_mj << "quadraple" << endl;
										done[i1]=true;done[i2]=true;done[i3]=true;done[i4]=true;
										npt[4]++;
										if(INPUT.system == "hydroxide")
										{
											if(last_nacc[i1+1]>=4) npt_hyper[4]++; // pay special attention to this
											else if(last_nacc[i1+1]==3) npt_tetra[4]++;
										}
										else if(INPUT.system == "hydronium")
										{
											if(last_nacc[i1] == 0) npt_hyper[4]++;
											else if(last_nacc[i1] == 1) npt_tetra[4]++;
										}
									}
								}
								else
								{
									cout << "possible quadraple: " << si[i1] << " total time " << total_time << endl;;
									quadraplej=true;
									ofs_mj << "quadraple" << endl;
									done[i1]=true;done[i2]=true;done[i3]=true;done[i4]=true;
									npt[4]++;
									if(INPUT.system == "hydroxide")
									{
										if(last_nacc[i1+1]>=4) npt_hyper[4]++; // pay special attention to this
										else if(last_nacc[i1+1]==3) npt_tetra[4]++;
									}
									else if(INPUT.system == "hydronium")
									{
										if(last_nacc[i1]==0) npt_hyper[4]++;
										else if(last_nacc[i1]==1) npt_tetra[4]++;
									}
								}
							}
						}// ai[i3]
					}// ai[i2]
				}// ai[i1]
			}// i6
		}

		// triple jump
		if(INPUT.func==1)
		{
			if(sixj==false and fivej==false and quadraplej==false and i4<ndim)
			{
				if(ai[i1] != ai[i3] and ai[i1] != ai[i4])
				{
					if(ai[i2] != ai[i4])
					{
						total_time = st[i4]-st[i2];
						if( total_time < dt_max )
						{
							cout << "triple jump: " << si[i1] << " total time " << total_time << endl;
							triplej=true;
							ofs_mj << "triple" << endl;
							done[i1]=true;done[i2]=true;done[i3]=true;
							npt[3]++;
							if(INPUT.system == "hydroxide")
							{
								if(last_nacc[i1+1]>=4) npt_hyper[3]++;
								else if(last_nacc[i1+1]==3) npt_tetra[3]++;
							}
							else if(INPUT.system == "hydronium")
							{
								if(last_nacc[i1]==0) npt_hyper[3]++;
								else if(last_nacc[i1]==1) npt_tetra[3]++;
							}
						}
					}
				}
			}
		}
		else if(INPUT.func==2)
		{
			if(sixj==false and fivej==false and quadraplej==false and i7<ndim)
			{
				if(ai[i1] != ai[i3] and ai[i1] != ai[i4])
				{
					if(ai[i2] != ai[i4])
					{
						total_time = st[i4]-st[i2];
						if( total_time < dt_max )
						{
							if(ai[i3]==ai[i5] and ai[i2]==ai[i6] and ai[i1]==ai[i7])
							{
								if( st[i7]-st[i4] > dt_max )
								{
									cout << "triple jump: " << si[i1] << " total time " << total_time << endl;
									triplej=true;
									ofs_mj << "triple" << endl;
									done[i1]=true;done[i2]=true;done[i3]=true;
									npt[3]++;
									if(INPUT.system == "hydroxide")
									{
										if(last_nacc[i1+1]>=4) npt_hyper[3]++;
										else if(last_nacc[i1+1]==3) npt_tetra[3]++;
									}
									else if(INPUT.system == "hydronium")
									{
										if(last_nacc[i1] == 0) npt_hyper[3]++;
										else if(last_nacc[i1] == 1) npt_tetra[3]++;
									}
									
								}
							}
							else
							{
								cout << "triple jump: " << si[i1] << " total time " << total_time << endl;
								triplej=true;
								ofs_mj << "triple" << endl;
								done[i1]=true;done[i2]=true;done[i3]=true;
								npt[3]++;
								if(INPUT.system == "hydroxide")
								{
									if(last_nacc[i1+1]>=4) npt_hyper[3]++;
									else if(last_nacc[i1+1]==3) npt_tetra[3]++;
								}
								else if(INPUT.system == "hydronium")
								{
									if(last_nacc[i1] == 0) npt_hyper[3]++;
									else if(last_nacc[i1] == 1) npt_tetra[3]++;
								}
							}
						}
					} // end ai[i2]
				}// end ai[i1]
			}//end triple
		}

		//---------------------
		// double jump
		//---------------------
		if(INPUT.func==1)
		{
			if(sixj==false and fivej==false and quadraplej==false and triplej==false and i3<ndim)
			{
				if(ai[i1] != ai[i3])
				{
					total_time = st[i3]-st[i2];
					if( total_time < dt_max )
					{
						cout << "double jump: " << si[i1] << " total time " << total_time << endl;
						doublej=true;
						ofs_mj << "double" << endl;
						done[i1]=true;done[i2]=true;
						npt[2]++;
						if(INPUT.system == "hydroxide")
						{
							if(last_nacc[i1+1]>=4) npt_hyper[2]++;
							else if(last_nacc[i1+1]==3) npt_tetra[2]++;
						}
						else if(INPUT.system == "hydronium")
						{
							if(last_nacc[i1] == 0) npt_hyper[2]++;
							else if(last_nacc[i1] == 1) npt_tetra[2]++;
						}
						
					}// end total time
				}// end ai
			}// end double jump
		}
		else if(INPUT.func==2)
		{
			// eliminate rattling in double jumps
			if(sixj==false and fivej==false and quadraplej==false and triplej==false and i5<ndim)
			{
				if(ai[i1] != ai[i3])
				{
					total_time = st[i3]-st[i2];
					if( total_time < dt_max )
					{
						if(ai[i2] == ai[i4] and ai[i1] == ai[i5])
						{
							if( st[i5]-st[i3] > dt_max )
							{
								cout << "double jump: " << si[i1] << " total time " << total_time << " a2==a4" << endl;
								doublej=true;
								ofs_mj << "double" << endl;
								done[i1]=true;done[i2]=true;
								npt[2]++;
								if(INPUT.system == "hydroxide")
								{
									if(last_nacc[i1+1]>=4) npt_hyper[2]++;
									else if(last_nacc[i1+1]==3) npt_tetra[2]++;
								}
								else if(INPUT.system == "hydronium")
								{
									if(last_nacc[i1] == 0) npt_hyper[2]++;
									else if(last_nacc[i1] == 1) npt_tetra[2]++;
								}
							}
						}
						else 
						{
							cout << "double jump: " << si[i1] << " total time " << total_time << " a2!=a4" << endl;
							doublej=true;
							ofs_mj << "double" << endl;
							done[i1]=true;done[i2]=true;
							npt[2]++;
							if(INPUT.system == "hydroxide")
							{
								if(last_nacc[i1+1]>=4) npt_hyper[2]++;
								else if(last_nacc[i1+1]==3) npt_tetra[2]++;
							}
							else if(INPUT.system == "hydronium")
							{
								if(last_nacc[i1] == 0) npt_hyper[2]++;
								else if(last_nacc[i1] == 1) npt_tetra[2]++;
							}
						}
					}// end total_time
				}// ai
			}// end 
		}// end func 2 

		// single
		if(sixj==false and fivej==false and quadraplej==false and triplej==false and doublej==false and i3<ndim)
		{
			if(ai[i1] != ai[i3])
			{
				cout << "single jump: " << si[i1] << endl;
				ofs_mj << "single" << endl;
				singlej = true;
				done[i1]=true;
				npt[1]++;
				if(INPUT.system == "hydroxide")
				{
					if(last_nacc[i1+1]>=4) npt_hyper[1]++;
					else if(last_nacc[i1+1]==3) npt_tetra[1]++;
				}
				else if(INPUT.system == "hydronium")
				{
					if(last_nacc[i1] == 0) npt_hyper[1]++;
					else if(last_nacc[i1] == 1) npt_tetra[1]++;
				}
			}
			else if(ai[i1] == ai[i3])
			{
				total_time = st[i3]-st[i2];
				if( total_time > dt_max)
				{
					cout << "single jump: " << si[i1] << endl;
					ofs_mj << "single" << endl;
					singlej = true;
					done[i1]=true;
					npt[1]++;
					if(INPUT.system == "hydroxide")
					{
						if(last_nacc[i1+1]>=4) npt_hyper[1]++;
						else if(last_nacc[i1+1]==3) npt_tetra[1]++;
					}
					else if(INPUT.system == "hydronium")
					{
						if(last_nacc[i1] == 0) npt_hyper[1]++;
						else if(last_nacc[i1] == 1) npt_tetra[1]++;
					}
					
				}
			}
		}

		if(sixj==false and fivej==false and quadraplej==false and triplej==false and doublej==false and singlej==false)
		{
			ofs_mj << "rattling" << endl;
			npt[0]++;
			if(INPUT.system == "hydroxide")
			{
				if(last_nacc[i1+1]>=4) npt_hyper[0]++;
				else if(last_nacc[i1+1]==3) npt_tetra[0]++;
			}
			else if(INPUT.system == "hydronium")
			{
				if(last_nacc[i1] == 0) npt_hyper[0]++;
				else if(last_nacc[i1] == 1) npt_tetra[0]++;
			}
		}


	}

	cout << setw(5) << "PT" << setw(10) << "number" << setw(15) << "via_hyper" << setw(15) << "via_tetra" << endl;
	ofs_running << setw(5) << "PT" << setw(10) << "number" << setw(15) << "via_hyper" << setw(15) << "via_tetra" << endl;
	for(int i=0; i<7; ++i)
	{
		cout << setw(5) << i << setw(10) << npt[i] << setw(15) << npt_hyper[i] << setw(15) << npt_tetra[i] << endl;
		ofs_running << setw(5) << i << setw(10) << npt[i] << setw(15) << npt_hyper[i] << setw(15) << npt_tetra[i] << endl;
	}


	cout << "For every 10 ps: " << endl;
	cout << setw(5) << "PT" << setw(10) << "number" << endl; 
	ofs_running << "For every 10 ps: " << endl;
	ofs_running << setw(5) << "JUMP" << setw(10) << "nPT" << setw(15) << endl; 
	double sum_pt=0.0;
	for(int i=0; i<7; ++i)
	{
		//cout << setw(5) << i << setw(10) << npt[i]/st[ndim-1]*10*i << endl;
		//ofs_running << setw(5) << i << setw(10) << npt[i]/st[ndim-1]*10*i << endl;
		cout << setw(5) << i << setw(10) << npt[i]/(st[ndim-1]-st[0])*10 << endl;
		ofs_running << setw(5) << i << setw(10) << npt[i]/(st[ndim-1]-st[0])*10 << endl;
		sum_pt+=i*npt[i]/(st[ndim-1]-st[0])*10;
	}
	cout << "Total number of PT events in every 10 ps: " << sum_pt << endl;
	ofs_running << "Total number of PT events in every 10 ps: " << sum_pt << endl;
		
	delete[] si;
	delete[] st;
	delete[] ai;
	delete[] nacc;
	delete[] ndon;

	delete[] last_si;
	delete[] last_st;
	delete[] last_ai;
	delete[] last_nacc;
	delete[] last_ndon;

	delete[] indexH;

	delete[] done;

	delete[] npt;
	delete[] npt_hyper;
	delete[] npt_tetra;

	ifs.close();
	ofs_mj.close();

	return;
}
