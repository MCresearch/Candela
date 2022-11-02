//------------------------------------------
// STRUCTURE OF CLASS:
//   CLASS SSF
//     |_FUNCTION Routine
//     |_FUNCTION check_file_exist
//     |_FUNCTION cal
//     |_FUNCTION ssf_3D
//       |_FUNCTION sumup
//     |_FUNCTION cal_diff_norm
//     |_FUNCTION ssf_1D
//     |_FUNCTION write_ssf
//------------------------------------------
#include "cellFile.h"
#include "input.h"
#include "ssf.h"
#include "math.h"

void SSF::Routine()
{
	TITLE("SSF","Routine");
	this->cal();

	return;
}

void SSF::check_file_exist(const string &name)
{
	// compute static structure factor
	ifstream ifs(name.c_str());
	if(ifs)
	{
		cout << " The SSF has been calculated, quit." << endl;
		exit(0);
	}
	else
	{
		cout << " The file " << name << " is calculating now!" << endl;
	}
}


//-------------------------------------------------------------
// we will calculate the static structure factor S(g) at each
// g point, where the interval of g points is delta_g,
//-------------------------------------------------------------
void SSF::cal()
{
	TITLE("SSF","cal");

	// --> INITIALIZE <--

	assert(INPUT.struf_dgx>0);
	assert(INPUT.struf_dgy>0);
	assert(INPUT.struf_dgz>0);
	assert(INPUT.struf_ng>0);

	// number of delta_g, determined by 2pi/a 
	// where a is the lattice constant.
	this->dgx = INPUT.struf_dgx;
	this->dgy = INPUT.struf_dgy;
	this->dgz = INPUT.struf_dgz;

	// number of g points along each direction.
	this->ngx = INPUT.struf_ng;
	this->ngy = INPUT.struf_ng;
	this->ngz = INPUT.struf_ng;

	// total nuber of g points (3 dimensional).
	this->ngtot = (2*ngx+1) * (2*ngy+1) * (2*ngz+1);

	// static structure factor in 3 dimensional.
	float* sf = new float[ngtot];
	ZEROS(sf, ngtot);
	cout << " Number of total K=2pi/L*(" 
		<< ngx << "," 
		<< ngy << "," 
		<< ngz << ") = " << ngtot <<endl;

	// ---> BODY <--- 

	// circle from geo_1 to geo_2.
	// because the final static structure factor
	// is an average property, so we need to have
	// a log of samples here.
	assert(INPUT.geo_2 >= INPUT.geo_1);
	assert(INPUT.geo_interval > 0);

	int count_geometry_number=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; igeo+=INPUT.geo_interval)
	{
		// cel_in : input geometry file.
		CellFile cel_in;

		stringstream ss; ss << igeo;
		cel_in.file_name = ss.str(); 
		// cout << " File name is " << ss.str() << endl;

		// Read in geometry.
		if( !CellFile::ReadGeometry( cel_in ) ) continue;
		if(igeo<INPUT.geo_ignore)//qinrui add 2020-1-6
        {
            cout<<"ignore:"<<igeo<<endl;
            cel_in.clean();
            continue;
        }
		++count_geometry_number;
#ifdef __MPI
		if((count_geometry_number-1)%NPROC != RANK)
        {
            cel_in.clean();
            continue;
        }
#endif
		cout<<"igeo="<<igeo<<endl;
		// (1) Calculate the static structure factor (3D).
		this->ssf_3D( cel_in, sf );
		cel_in.clean();//qianrui add 2020-1-6
	}
	
	// do average of 3D ssf.
	assert(count_geometry_number>0);
	cout << " count_geometry_number = " << count_geometry_number << endl;
	if(count_geometry_number>0)
	{
		for(int ig=0; ig<ngtot; ++ig)
		{
			sf[ig] /= count_geometry_number;
		}
	}

#ifdef __MPI
	MPI_Allreduce(MPI_IN_PLACE, sf, ngtot, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#endif
	
	if(RANK==0)
	{
		// (2) Prepare the G vectors. Next we will project 3D ssf
		// to 1D.
		int* nG_1D = new int[ngtot]; // Number of |G|.
		int* norm_index = new int[ngtot]; // Index between 3D G and 1D |G|.
		ZEROS( nG_1D, ngtot );
		ZEROS( norm_index, ngtot );
		this->diff_norm = cal_diff_norm( nG_1D, norm_index );

		// (3) Calculate the 1D SSF.
		float* G_1D = new float[diff_norm]; // Norm of |G|
		float* sf_1D = new float[diff_norm]; // SSF in 1D.
		ZEROS( G_1D, diff_norm );
		ZEROS( sf_1D, diff_norm );
		this->ssf_1D( G_1D, sf_1D, norm_index, nG_1D, sf );
		this->rank_ssf( diff_norm, G_1D, sf_1D, nG_1D);
		
		// (4) Output the final SSF.
		this->write_ssf( diff_norm, G_1D, sf_1D );
		this->write_smoothssf(G_1D, sf_1D, nG_1D, diff_norm, INPUT.struf_avdg, "sm-"+INPUT.ssf_out);

		delete[] nG_1D;
		delete[] norm_index;
		delete[] G_1D;
		delete[] sf_1D;
	}

	// --> CLEAN <--

	delete[] sf;
	return;
}

void SSF::ssf_3D( 
	const Cell &cel, // cell information
	float* sf // final structure factor
	) const
{
	TITLE("SSF","ssf_3D");

	// --> INITIALIZE <--
	//qianrui begin
	float* sum_cos = new float[ngtot];
	float* sum_sin = new float[ngtot];
    ZEROS( sum_cos, ngtot );
    ZEROS( sum_sin, ngtot );
	//qianrui end

    //// Calculate the distance between atoms.
    float x2, y2, z2; // atom positions for atom 2.
    float dr[3]; // delta x,y,z between atom1, atom2.

    float dis;
    int ito=-1;
    int ith=-1;
    int itc=-1;
    for(int it=0;it <INPUT.ntype; ++it)
    {
        if(cel.atom[it].id=="O") ito=it;
        else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
        else if(cel.atom[it].id=="C") itc=it;
    }

	// --> BODY <--
	//qianrui begin
    if(INPUT.func==1)
    {
        for(int it=0; it<INPUT.ntype; ++it)
        {
            for(int ia=0; ia<cel.atom[it].na; ++ia)
            {
	    				dr[0] = cel.atom[it].pos[ia].x;
	    				dr[1] = cel.atom[it].pos[ia].y;
	    				dr[2] = cel.atom[it].pos[ia].z;
    
    					this->sumup(sum_cos, sum_sin, dr);
            }
        }
    }
    else if(INPUT.func==2)
    {
	cout << INPUT.func << endl;
		if (INPUT.func_b == 1) 
		{
			for(int ia=0; ia<cel.atom[ito].na; ++ia)
			{
				dr[0] = cel.atom[ito].pos[ia].x;
				dr[1] = cel.atom[ito].pos[ia].y;
				dr[2] = cel.atom[ito].pos[ia].z;

				this->sumup(sum_cos, sum_sin, dr);
			}
		}
		if (INPUT.func_b == 2) // renxi added 20220329
		{
			int ik=0;
			for(int ix=-ngx; ix<=ngx; ++ix)
			{
				for(int iy=-ngy; iy<=ngy; ++iy)
				{
					for(int iz=-ngz; iz<=ngz; ++iz)
					{
						++ik;

					// because we only need to use half of the 
					// G vectors (same in gamma-only algorithm)
						if(ix<0) continue;
						if (ix == 0) sum_cos[ik] = 1;
						else sum_cos[ik] = 2;
					}
				}
			}
			for (int ia1=0; ia1<cel.atom[ito].na-1; ++ia1)
			{
				for (int ia2=ia1; ia2<cel.atom[ito].na; ++ia2)
				{
					dr[0] = cel.atom[ito].pos[ia1].x - cel.atom[ito].pos[ia2].x;
					dr[1] = cel.atom[ito].pos[ia1].y - cel.atom[ito].pos[ia2].y;
					dr[2] = cel.atom[ito].pos[ia1].z - cel.atom[ito].pos[ia2].z;
					this->sumup(sum_cos, sum_sin, dr);
				}
			}
		}
    }//renxi added 20201112
	//qianrui end

	// Use the final expression of structure factor.
	// s(G)=1 + 1/N * [ \sum_{ij,j!=i}exp(iGr) ] 
	if (INPUT.func_b == 1)
	{
		for(int ig=0; ig<ngtot; ++ig)
		{
			// 2.0 conunts for the atom pairs.
			double this_sf = 0.0;
			this_sf = pow(sum_cos[ig],2)+pow(sum_sin[ig],2);//qianrui;because there is only one index i, 2.0 counts doesn't exist. And because of it, sqrt(2) below is introduced.
			if(INPUT.func == 1)
			{
				this_sf/= INPUT.natom;
			}
			else if(INPUT.func == 2)
			{
				this_sf/= INPUT.natom1;
			}
			// mohan: be careful, 
			sf[ig] += this_sf;

			//cout << " ig=" << ig << " sum_exp=" << sum_exp[ig] << endl;
		}
	}
	if (INPUT.func_b == 2)
	{
		for(int ig=0; ig<ngtot; ++ig)
		{
			// 2.0 conunts for the atom pairs.
			double this_sf = 0.0;
			this_sf = sum_cos[ig];//qianrui;because there is only one index i, 2.0 counts doesn't exist. And because of it, sqrt(2) below is introduced.
			if(INPUT.func == 1)
			{
				this_sf/= INPUT.natom;
			}
			else if(INPUT.func == 2)
			{
				this_sf/= INPUT.natom1;
			}
			// mohan: be careful, 
			sf[ig] += this_sf;

			//cout << " ig=" << ig << " sum_exp=" << sum_exp[ig] << endl;
		}
	}
	// --> CLEAN <--

	delete[] sum_cos;
	delete[] sum_sin;

	return;
}


void SSF::sumup( float *sum_cos,float *sum_sin, const float dr[3] ) const
{
	int ik=0;
	if (INPUT.func_b == 1)
	{
		for(int ix=-ngx; ix<=ngx; ++ix)
		{
			for(int iy=-ngy; iy<=ngy; ++iy)
			{
				for(int iz=-ngz; iz<=ngz; ++iz)
				{
					++ik;

					// because we only need to use half of the 
					// G vectors (same in gamma-only algorithm)
					if(ix<0) continue;

					// For the diagonalize part, the factor is 1,
					// and for the non-diagonalize G vectors, it's 2.
					float factor=sqrt(2);//qianrui because there exists square power.
					if(ix==0) factor=1.0;

					// Get the G vector.
					float kx = ix * this->dgx;
					float ky = iy * this->dgy;
					float kz = iz * this->dgz;

					// exp(ik*(r_i - r_j)) appears here!
					float phase = kx * dr[0] + ky * dr[1] + kz * dr[2];
					sum_cos[ik-1] += factor*cos(phase);//qianrui
					sum_sin[ik-1] += factor*sin(phase);//qianrui
				}
			}
		}
	}
	if (INPUT.func_b == 2)
	{
		for(int ix=-ngx; ix<=ngx; ++ix)
		{
			for(int iy=-ngy; iy<=ngy; ++iy)
			{
				for(int iz=-ngz; iz<=ngz; ++iz)
				{
					++ik;
					if(ix<0) continue;
					float kx = ix * this->dgx;
					float ky = iy * this->dgy;
					float kz = iz * this->dgz;
					float phase = kx * dr[0] + ky * dr[1] + kz * dr[2];
					sum_cos[ik-1] += 4*cos(phase);
				}
			}
		}
	}// renxi 20220329
	return;
}

int SSF::cal_diff_norm(
	int *nG_1D, 
	int *norm_index
) const
{
	TITLE("SSF","cal_diff_norm");

	// --> INITIALIZE <--

	int* norm_tanker = new int[ngtot];
	ZEROS( norm_tanker, ngtot );
	int diff_norm = 0;
	for(int ig=0; ig<ngtot; ++ig)
	{
		norm_tanker[ig] = -1;
		nG_1D[ig] = 0;
		norm_index[ig] = -1;
	}

	// --> BODY <--

	int ig_global=0;
	for(int ix=-ngx; ix<=ngx; ++ix)
	{
		for(int iy=-ngy; iy<=ngy; ++iy)
		{
			for(int iz=-ngz; iz<=ngz; ++iz)
			{
				bool not_new = false;
				int norm = ix*ix+iy*iy+iz*iz;

				//// Compare to the old |G|.
				for(int i=0; i<diff_norm; ++i)
				{
					if( norm_tanker[i] == norm )
					{
						not_new = true;
						norm_index[ig_global] = i;
					}
				}

				if(not_new)
				{
					nG_1D[norm_index[ig_global]] += 1;
				}
				else
				{
					//// Find a new |G| !
					norm_tanker[diff_norm] = norm;
					nG_1D[diff_norm] += 1;
					norm_index[ig_global] = diff_norm;
					++diff_norm;
		//			cout << " The new norm is = " << sqrt((float)norm) << endl;
				}
				++ig_global;
			}
		}
	}
	cout << " Diff_norm = " << diff_norm<< endl;

	// --> CLEAN <--

	delete[] norm_tanker;
	return diff_norm;
}
	

void SSF::ssf_1D( 
	float *G_1D, 
	float *sf_1D, 
	const int *norm_index,
	const int *nG_1D,
	const float *sf	
) const
{
	TITLE("SSF","ssf_1D");

	// --> INITIALIZE <--

	for(int i=0; i<diff_norm; ++i)
	{
		G_1D[i] = -1.0;
		sf_1D[i] = 0.0;
	}

	// --> BODY <--

	int ig_global=0;
	for(int ix=-ngx; ix<=ngx; ++ix)
	{
		for(int iy=-ngy; iy<=ngy; ++iy)
		{
			for(int iz=-ngz; iz<=ngz; ++iz)
			{
				const int ig_1D = norm_index[ig_global];

                G_1D[ig_1D] = pow(ix*dgx,2)+pow(iy*dgy,2)+pow(iz*dgz,2);
                G_1D[ig_1D] = sqrt(G_1D[ig_1D]);

//				if(G_1D[ig_1D]<0.0)
//				{
//					G_1D[ig_1D] = sqrt((float)ix*ix+iy*iy+iz*iz)*dg;
//				}

				sf_1D[ig_1D] += sf[ig_global];
				++ig_global;
			}
		}	
	}


	for(int i=0; i<diff_norm; ++i)
	{
		assert(nG_1D[i]>0);
		sf_1D[i]/=nG_1D[i];
		//ofs << nor[i] << " " << nG_1D[i] << endl;
	}

	// --> CLEAN <--

	// no clean

	return;
}


void SSF::rank_ssf( const int diff_norm, float *G_1D, float* sf_1D, int* nG_1D ) const
{
	TITLE("SSF","rank_ssf");

	// --> INITIALIZE <--
	bool *visited = new bool[diff_norm];
	for(int i=0; i<diff_norm; ++i) visited[i] = false;
	float *tmp_G1D = new float [diff_norm];
	float *tmp_sf1D = new float [diff_norm];
	int *tmp_nG1D = new int [diff_norm];
	for(int i = 0 ; i < diff_norm ; ++i)
	{
		tmp_G1D[i] = G_1D[i];
		tmp_sf1D[i] = sf_1D[i];
		tmp_nG1D[i] = nG_1D[i];
	}


	// --> BODY <--
	for(int i=0; i<diff_norm; ++i)
	{
		int index = 0;
		float min = 10000000;
		for(int j=0; j<diff_norm; ++j)
		{
			if( visited[j] ) continue;
			else if( tmp_G1D[j] < min )
			{
				index = j;
				min = tmp_G1D[j];
			}
		}
		visited[index] = true;
		G_1D[i] = tmp_G1D[index];
		sf_1D[i] = tmp_sf1D[index];
		nG_1D[i] = tmp_nG1D[index];
	}

	// --> CLEAN <--

	delete[] visited;
	delete[] tmp_G1D;
	delete[] tmp_sf1D;
	delete[] tmp_nG1D;
	return;
}

void SSF::write_ssf( const int diff_norm, const float *G_1D, const float* sf_1D ) const
{
	assert(diff_norm > 0);
	// output the static structure factor.
	ofstream ofs(INPUT.ssf_out.c_str());
	for(int i = 0 ; i < diff_norm ; ++ i)
	{
		if(G_1D[i]>0)
			ofs << G_1D[i] << " " << sf_1D[i] << endl;
	}
}

void SSF::write_smoothssf(const float *G_1D,const float *sf, const int *nG_1D, const int diff_norm, const double dG, const string filename) const
{
	assert(dG > 0);
	ofstream ofs(filename.c_str());
	int i0 = 0;
	double tm_g = 0;
	double tm_ssf = 0;
	double tm_pre = 0;
	for(int i = 0 ; i < diff_norm ; ++i)
	{
		float G = G_1D[i];
		int nG = nG_1D[i];
		if(G < 1e-6) continue;
		if(int(G/dG) == i0)
		{
			tm_pre += nG;
			tm_g += G * nG;
			tm_ssf += sf[i] * nG; 
		}
		else
		{
			if(tm_g > 0)
			{
				tm_g /= tm_pre;
				tm_ssf /= tm_pre;
				ofs<<tm_g<<' '<<tm_ssf<<endl;
			}
			tm_pre = nG;
			tm_g = G * nG;
			tm_ssf = sf[i] * nG;
			i0 = int(G/dG); 
		}
	}
	//i == diff_norm
	if(tm_g > 0)
	{
		tm_g /= tm_pre;
		tm_ssf /= tm_pre;
		ofs<<tm_g<<' '<<tm_ssf<<endl;
	}

	ofs.close();
	return;
}
