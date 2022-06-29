//------------------------------------------
// STRUCTURE OF CLASS:
//   CLASS DSF
//     |_FUNCTION Routine
//     |_FUNCTION cal
//------------------------------------------
#include "cellFile.h"
#include "input.h"
#include "dsf.h"
#include "math.h"

void DSF::Routine()
{
	TITLE("DSF","Routine");
	this->cal();

	return;
}



void DSF::cal()
{
	TITLE("DSF","cal");

	// --> INITIALIZE <--
	//// g vectors information.
	assert(INPUT.struf_dgx>0);
	assert(INPUT.struf_ng>0);
	this->dg = INPUT.struf_dgx;
	cout << " not finished yet, please change dg" << endl;
	exit(0);
	this->ngx = INPUT.struf_ng;
	this->ngy = INPUT.struf_ng;
	this->ngz = INPUT.struf_ng;
	this->ngtot = (2*ngx+1) * (2*ngy+1) * (2*ngz+1);
	cout << " ngtot=" << ngtot << endl;
	assert(ngtot<10000);
	
	//// time information.
	assert(INPUT.geo_2 >= INPUT.geo_1);
	assert(INPUT.dsf_dt > 0);
	this->timetot = INPUT.dsf_dt * (INPUT.geo_2 - INPUT.geo_1 + 1);
	cout << " Total Time is " << timetot/1000 << " pico-second" << endl;
	this->nt = this->count_geo();	
	cout << " Number of differnt time: " << nt << endl;

	float* sum_exp = new float[ngtot];
	ZEROS(sum_exp, ngtot);
	
	float** sf_3D = new float*[nt];
	for(int it=0; it<nt; ++it)
	{
		sf_3D[it] = new float[ngtot];
		ZEROS(sf_3D[it], ngtot);
	}
	
	cout << " Number of total K=2pi/L*(" 
		<< ngx << "," 
		<< ngy << "," 
		<< ngz << ") = " << ngtot <<endl;

	// --> BODY <--- 

	const int nvel = INPUT.geo_2 - INPUT.geo_1 + 1;
	int count_neqi = 0;
	for(int ieqi=0; ieqi<nvel; ++ieqi)
	{
		bool file_exist = true;

		CellFile cel1;
		stringstream ss; ss << ieqi+INPUT.geo_1;
		cel1.file_name = ss.str();
		if(! CellFile::ReadGeometry( cel1 ))
		{
			file_exist = false;
			continue;
		}

		// if the number of equilibrium states
		// has reach the required number of states
		// 'neqi', then we stop including more 
		// equilibrium states.
		++count_neqi;
		if(count_neqi==INPUT.dsf_neqi) break;

		//// circle from geo_1 to geo_2.
		int time_now = 0;
		for(int igeo=0; igeo<nvel; ++igeo)
		{
			//// cel_in : input geometry file.
			CellFile cel2;

			stringstream ss; ss << igeo+INPUT.geo_1+ieqi;
			cel2.file_name = ss.str(); 
			//// cout << " File name is " << ss.str() << endl;

			//// Read in geometry.
			if( !CellFile::ReadGeometry( cel2 ) ) continue;

			////----------------------------------------------
			//// Calculate the static structure factor (3D).
			//// if cel1=cel2, the result should be 1,
			//// no matter what q it is.
			////----------------------------------------------
			this->dsf_3D( cel1, cel2, sum_exp, sf_3D[time_now] );

			++time_now;
		}
	}

	// (2)
	int* nG_1D = new int[ngtot];
	int* norm_index = new int[ngtot];
	this->diff_norm = cal_diff_norm( nG_1D, norm_index );

	// (3)
	float* G_1D = new float[diff_norm];
	float** sf_1D = new float*[nt];
	for(int it=0; it<nt; ++it)
	{
		sf_1D[it] = new float[diff_norm];
	}
	this->dsf_1D( G_1D, sf_1D, norm_index, nG_1D, sf_3D);

	// (4)
	this->write_dsf( G_1D, sf_1D );

	// --> CLEAN <---
	delete[] sum_exp;
	for(int it=0; it<nt; ++it)
	{
		delete[] sf_3D[it];
		delete[] sf_1D[it];
	}
	delete[] sf_3D;
	delete[] sf_1D;
	delete[] G_1D;
	delete[] norm_index;
	delete[] nG_1D;
	return;
}

void DSF::dsf_3D( 
	const Cell &cel1, // cell1 (time=t0)
	const Cell &cel2, // cell2 (time=t) 
	float* sum_exp, // calculate part of the structure factor
	float* sf_3D // final structure factor
	) const
{
	TITLE("DSF","ssf_3D");

	// --> INITIALIZE <--
	ZEROS(sum_exp, ngtot);

    //// Calculate the distance between atoms.
    float x2, y2, z2; // atom positions for atom 2.
    float dr[3]; // delta x,y,z between atom1, atom2.
    float dis;


	// --> BODY <--
    for(int it=0; it<INPUT.ntype; ++it)
    {
        for(int ia=0; ia<cel1.atom[it].na; ++ia)
        {
			//double check the atom positions.
			//cout << cel1.atom[it].pos[ia].x << " "
			//<< cel1.atom[it].pos[ia].y << " "
			//<< cel1.atom[it].pos[ia].z << endl;

            // Search in the around cells.
			// Try to find the shortest atom distance
			// between atom 1 and atom 2,
			// then that's the distance we want!
			float shortest_distance2 = 10000.0;
			int which_i, which_j, which_k;
			for(int i=-1; i<=1; ++i)
			{
				for(int j=-1; j<=1; ++j)
				{
					for(int k=-1; k<=1; ++k)
					{
						// add cell length to atom 2
						cel2.add_cell_length(it, ia, i, j, k, x2, y2, z2);
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
					}
				}
			}
			//cout << " " << which_i << " " << which_j << " " << which_k << " cell is for atom " << ia << endl; 
			// Here we identify the atom in cell: (which_i, which_j, which_k)
					
			// we get the vector 'dr' again.
			cel2.add_cell_length(it, ia, which_i, which_j, which_k, x2, y2, z2);
			dr[0] = cel1.atom[it].pos[ia].x - x2;
			dr[1] = cel1.atom[it].pos[ia].y - y2;
			dr[2] = cel1.atom[it].pos[ia].z - z2;

			this->sumup(sum_exp, dr);
		}
    }

	// Use the final expression of structure factor.
	// s(G)=1 + 1/N * [ \sum_{ij,j!=i}exp(iGr) ] 
	for(int ig=0; ig<ngtot; ++ig)
	{
		sf_3D[ig] = sum_exp[ig];
		sf_3D[ig]/= INPUT.natom;
//		cout << " sf_3D[" << ig << "]=" << sf_3D[ig] << endl;
	}

	return;
}


void DSF::sumup( float *sum_exp, const float dr[3] ) const
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

				// For the diagonalize part, the factor is 1,
				// and for the non-diagonalize G vectors, it's 2.
				float factor=2.0;
				if(ix==0) factor=1.0;

				// Get the G vector.
				float kx = ix * this->dg;
				float ky = iy * this->dg;
				float kz = iz * this->dg;

				// exp(ik*(r_i - r_j)) appears here!
				float phase = kx * dr[0] + ky * dr[1] + kz * dr[2];
				sum_exp[ik-1] += factor*cos(phase);
				//cout << " ik=" << ik << " sum_exp=" << sum_exp[ik-1] << endl;
			}
		}
	}
	
	return;
}

int DSF::cal_diff_norm(
	int *nG_1D, 
	int *norm_index
) const
{
	TITLE("DSF","cal_diff_norm");

	// --> INITIALIZE <--
	int* norm_tanker = new int[ngtot];
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
	

void DSF::dsf_1D( 
	float *G_1D, 
	float **sf_1D, 
	const int *norm_index,
	const int *nG_1D,
	float **sf	
) const
{
	TITLE("DSF","dsf_1D");

	// --> INITIALIZE <--
	for(int i=0; i<diff_norm; ++i)
	{
		G_1D[i] = -1.0;
	}

	for(int it=0; it<nt; ++it)
	{
		ZEROS(sf_1D[it], diff_norm);
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

//				if(G_1D[ig_1D]<0.0)
				{
					G_1D[ig_1D] = sqrt((float)ix*ix+iy*iy+iz*iz)*dg;
				}

				for(int it=0; it<nt; ++it)
				{
					sf_1D[it][ig_1D] += sf[it][ig_global];
				}
				++ig_global;
			}
		}	
	}

	for(int it=0; it<nt; ++it)
	{
		for(int ig=0; ig<diff_norm; ++ig)
		{
			assert(nG_1D[ig]>0);
			sf_1D[it][ig]/=nG_1D[ig];
			//ofs << nor[i] << " " << nG_1D[i] << endl;
		}
	}

	// --> CLEAN <--
	// no clean

	return;
}


void DSF::write_dsf( const float *G_1D, float** sf_1D ) const
{
	TITLE("DSF","write_ssf");

	// --> INITIALIZE <--
	assert(INPUT.dsf_neqi>0);
	assert(this->diff_norm > 0);
	// output the static structure factor.
	double unit_t = INPUT.dsf_dt/1000; // pico-second

	for(int ig=0; ig<diff_norm; ++ig)
	{
		// don't output g=0 term.
		if(G_1D[ig]>0)
		{
			stringstream ss;
			ss << "dsf_" << G_1D[ig] << ".txt";
			cout << " write dsf into file " << ss.str() << endl;
			ofstream ofs(ss.str().c_str());

			for(int it=0; it<nt; ++it)
			{
				ofs << it << " " << sf_1D[it][ig]/INPUT.dsf_neqi << endl;
			}
			ofs.close();
		}
	}

	// --> CLEAN <--
	return;
}




int DSF::count_geo()
{
	//// circle from geo_1 to geo_2.
	int count_geometry_number=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		//// cel_in : input geometry file.
		CellFile cel2;

		stringstream ss; ss << igeo;
		cel2.file_name = ss.str(); 
		//// cout << " File name is " << ss.str() << endl;

		//// Read in geometry.
		if( !CellFile::CheckGeometry( cel2 ) ) continue;
		++count_geometry_number;

	}
	cout << " count_geometry_number=" << count_geometry_number << endl;
	return count_geometry_number;
}



