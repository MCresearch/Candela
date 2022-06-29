//------------------------------------------
// STRUCTURE OF CLASS:
//   CLASS SSF_Selected
//     |_FUNCTION Routine
//     |_FUNCTION cal
//     |_FUNCTION ssf_3D
//       |_FUNCTION sumup
//     |_FUNCTION cal_diff_norm
//     |_FUNCTION ssf_1D
//     |_FUNCTION write_ssf
//------------------------------------------
#include "cellFile.h"
#include "input.h"
#include "ssf_selected.h"
#include "math.h"

void SSF_Selected::Routine()
{
	TITLE("SSF_Selected","Routine");

	ifstream ifs2("SSF.input1");
	if(ifs2)
	{
		this->select_g(ifs2);
		return;
	}

	ifstream ifs("SSF.input");
	if(!ifs)
	{
		this->generate_input();
	}
	else
	{
		ofs_running << "The input file SSF.input exists, let's calculte SSFs." << endl; 
		this->cal(ifs);
	}

	return;
}

void SSF_Selected::generate_input()
{
	TITLE("SSF_Selected","generate_input");

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

	ofstream ofs("SSF.input0");


	float* norm0 = new float[ngtot]; 
	for(int i=0; i<ngtot; ++i)
	{
		norm0[i]=-1.0;
	}
	int countDiffNorm=0;

    for(int ix=-ngx; ix<=ngx; ++ix)
    {
		if(ix<0) continue;
		float x = this->dgx * ix;
        for(int iy=-ngy; iy<=ngy; ++iy)
        {
			float y = this->dgy * iy;
            for(int iz=-ngz; iz<=ngz; ++iz)
            {
				float z = this->dgz * iz;
				float norm = sqrt(x*x+y*y+z*z);

				bool find_same_norm=false;
				for(int i=0; i<countDiffNorm; ++i)
				{
					if(norm0[i] == norm) 
					{
						find_same_norm=true;
						break;
					}
				}
				if(find_same_norm==false)
				{
					norm0[countDiffNorm]=norm;
					countDiffNorm++;
					ofs << ix*dgx << " " << iy*dgy 
					<< " " << iz*dgz << " " << norm << endl;	
				}
			}
		}
	}
	
	ofs.close();

	delete[] norm0;

	cout << "SORT THE 4th COLUMN IN SSF.input0 to SSF.input1 AND RUN D310 AGAIN TO GET REDUCED G VECTORS." << endl;
	cout << "REMEMBER TO PUT THE ROW NUMBER ON TOP OF THE SSF.input1 FILE." << endl;

	return;
}

void SSF_Selected::select_g(ifstream &ifs2)
{
	int ng1d=0;
	ifs2 >> ng1d;
	cout << "The number of g points is " << ng1d << endl;

	ofstream ofs2("SSF.input");
	float x,y,z,norm;

	float norm0 = 0.025;
	float norm00 = 3.00;
	float norm000 = 8.00;
	float dnorm1 = 0.05;
	float dnorm2 = 0.10;
	float dnorm3 = 0.25;
	for(int i=0; i<ng1d; ++i)
	{
		ifs2 >> x >> y >> z >> norm;
		if(norm <= 3.0)
		{
			if(norm > norm0)
			{
				ofs2 << setw(15) << x << setw(15) << y << setw(15) << z << setw(25) << norm << endl;	
				norm0+=dnorm1;
			}
		}
		else if(norm > 3.0 && norm <= 8.00)
		{
			if(norm > norm00)
			{
				ofs2 << setw(15) << x << setw(15) << y << setw(15) << z << setw(25) << norm << endl;	
				norm00+=dnorm2;
			}
		}
		else if(norm > 8.00)
		{
			if(norm > norm000)
			{
				ofs2 << setw(15) << x << setw(15) << y << setw(15) << z << setw(25) << norm << endl;	
				norm000+=dnorm3;
			}
		}
		
	}
	
	ofs2.close();

	cout << "RUN SSF.input1 TO OBTAIN SSF.input, REMEMBER TO ADD NUMBER OF G POINTS IN SSF.input." << endl;
	cout << "THEN DELETE BOTH SSF.input0 and SSF.input1 AND RUN D310 AGAIN." << endl;

	return;
}


void SSF_Selected::cal(ifstream &ifs)
{
	ifs >> this->ngtot;
	ofs_running << " number of G points is " << ngtot << endl;

	gx = new float[ngtot];
	gy = new float[ngtot];
	gz = new float[ngtot];
	norm_value = new float[ngtot];
	for(int ig=0; ig<ngtot; ++ig)
	{
		ifs >> gx[ig] >> gy[ig] >> gz[ig] >> norm_value[ig];		
	}

	// static structure factor in 3 dimensional.
	float* sf = new float[ngtot];
	ZEROS(sf, ngtot);

	assert(INPUT.geo_2 >= INPUT.geo_1);
	assert(INPUT.geo_interval > 0);


	int count_geometry_number=0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; igeo+=INPUT.geo_interval)
	{
		// cel_in : input geometry file.
		CellFile cel_in;

		stringstream ss; ss << igeo;
		cel_in.file_name = ss.str();
		 cout << " File name is " << ss.str() << endl;
	
		// Read in geometry.
		if( !CellFile::ReadGeometry( cel_in ) ) continue;
		++count_geometry_number;

		// (1) Calculate the static structure factor (3D).
		this->ssf_3D( cel_in, sf );
	}

    // do average of 3D ssf.
    assert(count_geometry_number>0);
    ofs_running << " count_geometry_number = " << count_geometry_number << endl;
    if(count_geometry_number>0)
    {
        for(int ig=0; ig<ngtot; ++ig)
        {
            sf[ig] /= count_geometry_number;
        }
    }

#ifdef __MPI
    float* sf_loc = new float[ngtot];
    for(int ig=0; ig<ngtot; ++ig)
    {
        //ofs_running << ig << " " << sf[ig] << endl;
        sf_loc[ig] = sf[ig];
        sf[ig] = 0.0;
    }
    MPI_Allreduce(sf_loc, sf, ngtot, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    delete[] sf_loc;
#endif

	if(RANK==0)
	{
		ofstream ofs(INPUT.ssf_out.c_str());
		for(int ig=0; ig<ngtot; ++ig)
		{
			ofs << setw(25) << norm_value[ig] << setw(25) << sf[ig] << endl;
		}	
		ofs.close();
	}	

	delete[] gx;
	delete[] gy;
	delete[] gz;
	delete[] norm_value;
	delete[] sf;

	return;
}

void SSF_Selected::ssf_3D(
    const Cell &cel, // cell information
    float* sf // final structure factor
    ) const
{
    TITLE("SSF_Selected","ssf_3D");

    // --> INITIALIZE <--

    float* sum_exp = new float[ngtot];
    ZEROS( sum_exp, ngtot );

    //// Calculate the distance between atoms.
    float x2, y2, z2; // atom positions for atom 2.
    float dr[3]; // delta x,y,z between atom1, atom2.
    float dis;
	
	// --> BODY <--
    for(int it=0; it<INPUT.ntype; ++it)
    {
		cout << "natom is " << cel.atom[it].na << endl;
        for(int ia=0; ia<cel.atom[it].na; ++ia)
        {
            //// Double check the atom positions.
            //// cout << cel.atom[it].pos[ia].x << " "
            //// << cel.atom[it].pos[ia].y << " "
            //// << cel.atom[it].pos[ia].z << endl;

            //// it2 start from species 1, because its about 'atom pairs',
            //// so we only need to calculate once.
            int count_ia2=0;
            for(int it2=it; it2<INPUT.ntype; ++it2)
            //for(int it2=0; it2<INPUT.ntype; ++it2)
            {
                for(int ia2=ia; ia2<cel.atom[it].na; ++ia2)
                //for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
                {
                    count_ia2++;

                    // In ssf, for the identical atoms,
                    // they will contribute 1 in the final
                    // expression of ssf.
                    if(it==it2 && ia==ia2) continue;

#ifdef __MPI
                    if( (count_ia2-1) % NPROC != RANK) continue;
#endif

                    dr[0] = cel.atom[it].pos[ia].x - cel.atom[it2].pos[ia2].x;
                    dr[1] = cel.atom[it].pos[ia].y - cel.atom[it2].pos[ia2].y;
                    dr[2] = cel.atom[it].pos[ia].z - cel.atom[it2].pos[ia2].z;

                    this->sumup(sum_exp, dr);
                }
            }
            // For check,
            //cout << " ia" << ia << " count=" << count_ia2 << endl;
        }
    }

    // Use the final expression of structure factor.
    // s(G)=1 + 1/N * [ \sum_{ij,j!=i}exp(iGr) ]
    for(int ig=0; ig<ngtot; ++ig)
    {
        // 2.0 conunts for the atom pairs.
        double this_sf = 0.0;
        this_sf = 2.0*sum_exp[ig];
        this_sf/= INPUT.natom;
        // mohan: be careful,
#ifdef __MPI
        this_sf+= 1.0/(double)NPROC;
#else
        this_sf+= 1.0;
#endif
        sf[ig] += this_sf;

        //cout << " ig=" << ig << " sum_exp=" << sum_exp[ig] << endl;
    }

    // --> CLEAN <--

    delete[] sum_exp;

    return;
}

void SSF_Selected::sumup( float *sum_exp, const float dr[3] ) const
{
	for(int ik=0; ik<this->ngtot; ++ik)
	{
		// exp(ik*(r_i - r_j)) appears here!
		float phase = gx[ik] * dr[0] + gy[ik] * dr[1] + gz[ik] * dr[2];
		sum_exp[ik] += cos(phase);
	}
    return;
}
	
