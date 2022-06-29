//------------------------------------------
// STRUCTURE OF CLASS:
//   CLASS SSF
//     |_FUNCTION Routine
//     |_FUNCTION selectq
//     |_FUNCTION cal_Fqt
//       |_FUNCTION cal_eiqR_t
//       |_FUNCTION cal_eiqR_0
//     |_FUNCTION output_Fqt
//------------------------------------------
#include "cellFile.h"
#include "input.h"
#include "isf.h"
#include "math.h"

ISF::ISF()
{
	qsaved_x = new float[1];
	qsaved_y = new float[1];
	qsaved_z = new float[1];
}

ISF::~ISF()
{
	delete[] qsaved_x;
	delete[] qsaved_y;
	delete[] qsaved_z;
}

void ISF::Routine()
{
	TITLE("ISF","Routine");
	cout << " Calculate the intermediate scattering function " << endl;

	//--------------------------------------------------
	// Intermediate scattering function 
	// F(q,t) = 1/N * <
	// ( \sum_{j=1,N} e^{-iqR_j(t)}) *
	// ( \sum_{i=1,N} e^{ iqR_i{0}})
	// >
	//--------------------------------------------------

	// (float) q:
	this->target_q = INPUT.isf_target_q;
	// (int) <>:
	this->nconfig = INPUT.isf_nconfig;
	// (int) t:
	this->ncorrelation = INPUT.isf_ncorrelation;
	// (int) N:
	this->nat = INPUT.natom;
	// (int) select atom number for i
	this->nat1 = INPUT.natom1;
	// (int) select atom number for j
	this->nat2 = INPUT.natom2;
	// (string) output file name
	string out_file = INPUT.isf_outfile;

	assert(nconfig>0);
	assert(ncorrelation>0);
	assert(nat>0);
	assert(nat1>0);

	// F(q,t)
	float* Fqt = new float[ncorrelation];
	for(int i=0; i<ncorrelation; ++i) Fqt[i]=0.0;

	this->selectq();

	this->cal_Fqt(Fqt,qsaved_x,qsaved_y,qsaved_z);
	this->output_Fqt(out_file,Fqt);


	delete[] Fqt;

	return;
}

void ISF::selectq()
{
	TITLE("ISF","selectq");

	// -->> INITIALIZE <<--

	// number of delta_g, determined by 2pi/a
	// where a is the lattice constant.
	float dg = INPUT.isf_dg;
	cout << " dg=" << dg << endl; 

	// number of g points along each direction.
	int ngx = INPUT.isf_ng;
	int ngy = INPUT.isf_ng;
	int ngz = INPUT.isf_ng;

	assert(INPUT.isf_ng > 0);
	assert(INPUT.isf_dg > 0);

	// -->> FUNCTION <<--
	
	// how many q vectors that are close to
	// target_q
	this->count_q = 0;
	for(int ix=-ngx; ix<=ngx; ++ix)
	{
		for(int iy=-ngy; iy<=ngy; ++iy)
		{
			for(int iz=-ngz; iz<=ngz; ++iz)
			{
				bool not_new = false;
				int norm = ix*ix+iy*iy+iz*iz;
				//// Compare to the old |G|.
							
				float q = sqrt(norm)*dg; 
				if( abs(q- this->target_q) < 0.01 ) 
				{
					++count_q;
			//		cout << " q=" << q << endl;
				}
			}
		}
	}
	
	cout << " Target_q " << target_q << endl;
	cout << " Find " << count_q << " q vectors in G space." << endl;
	assert(count_q>0);

	delete[] qsaved_x;
	delete[] qsaved_y;
	delete[] qsaved_z;
	this->qsaved_x = new float[count_q];
	this->qsaved_y = new float[count_q];
	this->qsaved_z = new float[count_q];
	for(int iq=0; iq<count_q; ++iq)
	{
		this->qsaved_x[iq] = 0.0;
		this->qsaved_y[iq] = 0.0;
		this->qsaved_z[iq] = 0.0;
	}

	this->count_q = 0;
	for(int ix=-ngx; ix<=ngx; ++ix)
	{
		for(int iy=-ngy; iy<=ngy; ++iy)
		{
			for(int iz=-ngz; iz<=ngz; ++iz)
			{
				bool not_new = false;
				int norm = ix*ix+iy*iy+iz*iz;
				//// Compare to the old |G|.
							
				float q = sqrt(norm)*dg; 
				if( abs(q- this->target_q) < 0.01 ) 
				{
					qsaved_x[count_q] = ix * dg;
					qsaved_y[count_q] = iy * dg;
					qsaved_z[count_q] = iz * dg;
					++count_q;
			//		cout << " q=" << q << endl;
				}
			}
		}
	}
	
	return;
}


void ISF::cal_Fqt(float *Fqt, float *qx, float *qy, float *qz)
{
	TITLE("ISF","cal_Fqt");

	assert(INPUT.step_interval_dynamics > 0);
	assert(INPUT.isf_config_start > 0);

	// for atom1, get R_I(t0) in cel1.
	for(int is=1; is<=nconfig; ++is)
	{
		// we choose step_interval_dynamics in order to avoid
		// correlation between configurations that are close. 
		const int file1 = is * INPUT.step_interval_dynamics + INPUT.isf_config_start - 1;
		CellFile cel1;

		cout << "\n Ensemble " << is+1 << "/" << nconfig << endl;	
		cout << " should read in file ion." << file1 << ".dat" << endl;
		//cout << "should read in file ion." << file1 << endl; 

		stringstream ss;
		ss << file1;
		cel1.file_name=ss.str();

		if(!CellFile::ReadGeometry(cel1)) continue;

		// for atom2, get R_J(t+t0) in cel2.
		int ic=0;
		int count_file=0;
		while(ic<ncorrelation)
		{
			const int file2 = count_file + file1;
			CellFile cel2;

			stringstream ss2;
			ss2 << file2;
			cel2.file_name=ss2.str();

			if(!CellFile::ReadGeometry(cel2)) 
			{
				++count_file;
				// too many files are missing!
				if(count_file>100000)
				{
					break;
				}
				continue;
			}

			cout << " is=" << is+1 << " ic=" << ic+1 << " file1=" << file1 << " file2=" << file2 << endl;

			Fqt[ic] += cal_eiqR(qx,qy,qz,cel1,cel2); 

			++ic;
			count_file += INPUT.isf_dcorrelation;
		}
	}

	// average over configurations
	for(int ic=0; ic<ncorrelation; ++ic)
	{
		Fqt[ic]/=nconfig;
	}
	return;
}


float ISF::cal_eiqR(float *qx, float *qy, float *qz,
	const Cell &cel1, const Cell &cel2)
{
	float x2, y2, z2;// atom positions for atom 2.
	float dr[3]; // delta x,y,z between atom1, atom2.
	float dis;
	float sum_exp = 0.0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		if(cel1.atom[it].id != INPUT.ele1) continue;

		for(int ia=0; ia<cel1.atom[it].na; ++ia)
		{
			//// Double check the atom positions.
			//// cout << cel.atom[it].pos[ia].x << " "
			//// << cel.atom[it].pos[ia].y << " "
			//// << cel.atom[it].pos[ia].z << endl;
			int count2 = 0;

			//// it2 start from species 1, because its about 'atom pairs',
			//// so we only need to calculate once.
			for(int it2=0; it2<INPUT.ntype; ++it2)
			{
				if(cel1.atom[it2].id != INPUT.ele2) continue;

				for(int ia2=0; ia2<cel2.atom[it2].na; ++ia2)
				{
					float shortest_distance2 = 10000.0;
					int which_i, which_j, which_k;
					for(int i=-1; i<=1; ++i)
					{
						for(int j=-1; j<=1; ++j)
						{
							for(int k=-1; k<=1; ++k)
							{
								// add cell length to atom 2
								cel2.add_cell_length(it2, ia2, i, j, k, x2, y2, z2);
								// calculate the distance between two atoms |r_1 - r_2|
								dr[0] = cel1.atom[it].pos[ia].x - x2;
								dr[1] = cel1.atom[it].pos[ia].y - y2;
								dr[2] = cel1.atom[it].pos[ia].z - z2;
								// to save the calculation, we avoid using sqrt.
								dis = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];
								if(dis < shortest_distance2)
								{
									shortest_distance2=dis;
									which_i=i;
									which_j=j;
									which_k=k;
								}
							}//end k
						}//end j
					}//end i
					// Here we identify the atom in cell: (which_i, which_j, which_k)
					// we get the vector 'dr' again.
					cel2.add_cell_length(it2, ia2, which_i, which_j, which_k, x2, y2, z2);
					dr[0] = cel1.atom[it].pos[ia].x - x2;
					dr[1] = cel1.atom[it].pos[ia].y - y2;
					dr[2] = cel1.atom[it].pos[ia].z - z2;

					// exp ( -iqR_t ) * exp ( iqR_0)
					// = exp ( -iq (R_t-R_0) )
					for(int iq=0; iq<this->count_q; ++iq)
					{
						float phase = qx[iq] * dr[0] + qy[iq] * dr[1] + qz[iq] * dr[2];
						sum_exp += cos( phase );
					}

					++count2;
				}// end ia2
			}// end it2
		}// end ia1
	}// end it1

	// mohan add 2014-06-01
	if(INPUT.ntype==2)
	{
		// !!!!!!!!!!!!!!!!!!
		// be careful of cel1 and cel2! mohan note 2016-11-09
		int n1=0;
		int n2=0;
		for(int it=0; it<INPUT.ntype; ++it)
		{
			if(cel1.atom[it].id == INPUT.ele1) n1=cel1.atom[it].na;
			if(cel2.atom[it].id == INPUT.ele2) n2=cel2.atom[it].na;
		}
		assert(n1>0);
		assert(n2>0);
		sum_exp *= (float)INPUT.natom / n1 / n2; 
	}
	else if(INPUT.ntype==1)
	{
		sum_exp /= (float)INPUT.natom;
	}

	sum_exp /= this->count_q;


//	cout << " sum_exp = " << sum_exp << endl;
	return sum_exp;
}






// output the Fqt array.
void ISF::output_Fqt(const string &out_file, const float *Fqt)const
{
	ofstream ofs(out_file.c_str());
	if(!ofs)
	{
		cout << " Can not open file : " << out_file << endl;
		exit(0);
	}

	for(int i=0; i<this->ncorrelation; ++i)
	{
		ofs << i+1 << " " << Fqt[i] << endl;
	}

	ofs.close();
}

