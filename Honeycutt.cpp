#include "cellFile.h"
#include "input.h"
#include "Honeycutt.h"
#include "math.h"

Honeycutt::Honeycutt()
{
	this->max_num_adj = 100;
	this->max_name = 100;
}

Honeycutt::~Honeycutt()
{

}


// general subroutine for HoneycuttAnderson Algorithm.
void Honeycutt::Routine()
{
	TITLE("Honeycutt","Routine");

	cout << " Cal Honeycutt Anderson Short-Range Order." << endl;

	cal();

	return;
}


// general subroutine.
void Honeycutt::cal()
{
	TITLE("Honeycutt","cal");

	assert(INPUT.natom>0);

	// adj_index
	this->nadj = new int[INPUT.natom]();
	this->adj_index = new int*[INPUT.natom];
	for(int i=0; i<INPUT.natom; ++i)
	{
		adj_index[i] = new int[max_num_adj];
		for(int j=0; j<max_num_adj; ++j)
		{
			adj_index[i][j] = -1;
		}	
	}

	// HA name
	this->HA_name = new string[max_name]; 
	this->HA_percent = new double[max_name];
	for(int i=0; i<this->max_name; ++i)
	{
		HA_name[i]="none";
		HA_percent[i] = 0.0;
	}	

	// if bdf_movie==true
	if(INPUT.bdf_movie==1)
	{
		ofs_snapshot.open("snapshot.xyz");	
		assert(INPUT.natom>0);
		selected_atom = new bool[INPUT.natom];
	}
	
	assert(INPUT.geo_interval>0);
    int count_geometry_number=0;
	
	cout << " geo_1 is " << INPUT.geo_1 << endl;
	cout << " geo_2 is " << INPUT.geo_2 << endl;

    for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; igeo=igeo+INPUT.geo_interval)
    {
//		if(igeo%INPUT.geo_interval!=0) continue;

        CellFile cel;

		// get the file name
        stringstream ss; ss << igeo;
        cel.file_name = ss.str();

        // cel : read in geometry file
        if( !CellFile::ReadGeometry( cel ) ) continue;
        ++count_geometry_number;

		HA_pairs(cel);
	}

	// statics
	
	double sum=0.0;
	for(int i=0; i<max_name; ++i)
	{
		sum += HA_percent[i];
	}
	cout << "sum = " << sum << endl;
	assert(sum>0.0);
	

	ofs_running << setw(10) << "HA_name" << setw(15) << "number" << setw(15) << "percent(%)" << endl;
	for(int i=0; i<max_name; ++i)
	{
		if(HA_name[i]!="none")
		{
			ofs_running << setw(10) << HA_name[i] << setw(15) << HA_percent[i] << setw(15) << HA_percent[i]/sum*100;
			if(HA_name[i]=="155") ofs_running << setw(15) << "ico";
			else if(HA_name[i]=="154") ofs_running << setw(15) << "distorted-ico";
			else if(HA_name[i]=="143") ofs_running << setw(15) << "distorted-ico";
			else if(HA_name[i]=="144") ofs_running << setw(15) << "bcc";
			else if(HA_name[i]=="166") ofs_running << setw(15) << "bcc";
			else if(HA_name[i]=="142") ofs_running << setw(15) << "fcc or hcp";

			ofs_running << endl;
		}
	}


	// delete adj_index
	delete[] nadj;
	for(int i=0; i<INPUT.natom; ++i)
	{
		delete[] adj_index[i];
	}
	delete[] adj_index;
	delete[] HA_name;

	if(INPUT.bdf_movie==1)
	{
		ofs_snapshot.close();
		delete[] selected_atom;
	}

	return;
}

void Honeycutt::HA_pairs(const Cell &cel)
{
	TITLE("Honeycutt","HA_pairs");

//	cout << "Cell" << endl;
//	cout << cel.a1.norm() << " " << cel.a2.norm() << " " << cel.a3.norm() << endl;

	Honeycutt::setup_nadj(cel,nadj,adj_index,max_num_adj);

	if(INPUT.bdf_movie==1)
	{
		for(int iat=0; iat<INPUT.natom; ++iat)
		{
			selected_atom[iat]=false;
		}
	}

	// determine the parameters in HA definition
	for(int iat=0; iat<INPUT.natom; ++iat)
	{
		for(int in=0; in<this->nadj[iat]; ++in)
		{
			const int iat2 = this->adj_index[iat][in];
			shared_neighbours(iat,iat2,nadj,adj_index,max_num_adj);
		}
	}	

	if(INPUT.bdf_movie==1)
	{
		int ntotal=0;
		for(int iat=0; iat<INPUT.natom; ++iat)
		{
			if(selected_atom[iat]==true)
			{
				++ntotal;
			}
		}
		ofs_snapshot << ntotal << endl;
		ofs_snapshot << "snpashot of HA" << endl;

		int iat0=0;
		for(int it=0; it<INPUT.ntype; ++it)
		{
			for(int ia=0; ia<cel.atom[it].na; ++ia)
			{
				if(selected_atom[iat0]==true)
				{
					ofs_snapshot << cel.atom[it].id << " " << cel.atom[it].pos[ia].x 
					<< " " << cel.atom[it].pos[ia].y << " " << cel.atom[it].pos[ia].z << endl;
				}
				++iat0;
			}
		}
	}

	return;
}

void Honeycutt::setup_nadj(const Cell &cel, int* nadj_in, int** adj_index_in, const int &max_num_adj_in)
{
	for(int i=0; i<INPUT.natom; ++i){ nadj_in[i]=0; }
	for(int i=0; i<INPUT.natom; ++i)
	{
		for(int j=0; j<max_num_adj_in; ++j)
		{
			adj_index_in[i][j]=-1;
		}
	}

	// search for neighbors
	int iat=0;
	for(int it=0; it<INPUT.ntype; ++it)
	{
		for(int ia=0; ia<cel.atom[it].na; ++ia)
		{
			int iat2=0;
			for(int it2=0; it2<INPUT.ntype; ++it2)
			{
				for(int ia2=0; ia2<cel.atom[it2].na; ++ia2)
				{
					if(iat==iat2) {++iat2; continue;}
					double dis = distance(cel.atom[it].pos[ia],cel.atom[it2].pos[ia2],
							cel.a1.norm(),cel.a2.norm(),cel.a3.norm());
					if(dis<INPUT.rcut)
					{
						const int index = nadj_in[iat];
						adj_index_in[iat][index]=iat2;
						++nadj_in[iat];
						assert(nadj_in[iat]<max_num_adj_in);
					}
					++iat2;
				}
			}
//			cout << INPUT.natom << " " << iat2 << endl;
			assert(INPUT.natom==iat2);
//			cout << setw(10) << "Atom" << setw(10) << "Nadj" << endl;
//			cout << setw(10) << iat  << setw(10) << nadj_in[iat] << endl;
			for(int ii=0; ii<nadj_in[iat]; ++ii)
			{
//				cout << setw(10) << ii+1 << setw(10) << adj_index[iat][ii] << endl;
			}
 
			++iat;
		}
	}
	assert(INPUT.natom==iat);

	return;
}


// nsn: number of shared neighbours
// nsb: number of shared bonds
bool Honeycutt::search_neighbours(const int &iat, const int &iat2, int* nadj_in, int** adj_index_in, const int &max_num_adj_in,
const int &nsn, const int &nsb)
{
	assert(nsn<7);
	assert(nsb<7);

	if(iat==iat2) return false;

	bool is_neighbour=false;
	for(int in=0; in<nadj_in[iat]; ++in)
	{
		const int iat2_tmp = adj_index_in[iat][in];
		if(iat2==iat2_tmp)
		{
			is_neighbour=true;
		}
	}

	if(is_neighbour==false) return false;

	// list of shared atoms
	int* nnn = new int[max_num_adj_in];
	for(int i=0; i<max_num_adj_in; ++i)
	{
		nnn[i]=-1;
	}

	// second number
	int n_comm_neighbour=0;
	for(int ii=0; ii < nadj_in[iat2]; ++ii)
	{
		const int iat3 = adj_index_in[iat2][ii];
		for(int jj=0; jj < nadj_in[iat]; ++jj)
		{
			const int iat4 = adj_index_in[iat][jj];
			if(iat3==iat4)
			{
				nnn[n_comm_neighbour]=iat3;
				n_comm_neighbour++;
				//cout << iat3 << endl;
			}
		}
	} 

	if(n_comm_neighbour!=nsn) 
	{
		delete[] nnn;
		return false;
	}

	// third number
	int nbonds=0;
	for(int ii=0; ii<n_comm_neighbour; ++ii)
	{
		const int iat3 = nnn[ii];
		for(int kk=ii+1; kk<n_comm_neighbour; ++kk)
		{
			const int iat5 = nnn[kk];
			for(int jj=0; jj < nadj_in[iat3]; ++jj)
			{
				const int iat6 = adj_index_in[iat3][jj];
				if(iat5 == iat6)
				{
					++nbonds;
					break;
				}
			}
		}
	}
	
	if(nbonds!=nsb) 
	{
		delete[] nnn;
		return false;
	}

	delete[] nnn;

	return true;
}


// search for neighbours and do statistics
void Honeycutt::shared_neighbours(const int &iat, const int &iat2, int* nadj_in, int** adj_index_in, const int &max_num_adj_in)
{
	assert(iat!=iat2);

	// list of shared atoms
	int* nnn = new int[max_num_adj_in];
	for(int i=0; i<max_num_adj_in; ++i)
	{
		nnn[i]=-1;
	}


	// second number
	int n_comm_neighbour=0;
	for(int ii=0; ii < nadj_in[iat2]; ++ii)
	{
		const int iat3 = adj_index_in[iat2][ii];
		for(int jj=0; jj < nadj_in[iat]; ++jj)
		{
			const int iat4 = adj_index_in[iat][jj];
			if(iat3==iat4)
			{
				nnn[n_comm_neighbour]=iat3;
				n_comm_neighbour++;
				//cout << iat3 << endl;
			}
		}
	} 
	if(n_comm_neighbour>0)
	{
		//cout << setw(6) << "pair" << setw(8) << iat << setw(8) << iat2 << setw(10) << n_comm_neighbour << endl;
	}


	if(n_comm_neighbour<7)
	{
		// third number
		int nbonds=0;
		for(int ii=0; ii<n_comm_neighbour; ++ii)
		{
			const int iat3 = nnn[ii];
			for(int kk=ii+1; kk<n_comm_neighbour; ++kk)
			{
				const int iat5 = nnn[kk];
				for(int jj=0; jj < nadj_in[iat3]; ++jj)
				{
					const int iat6 = adj_index_in[iat3][jj];
					if(iat5 == iat6)
					{
						++nbonds;
						break;
					}
				}
			}
		}


		// please DIY to get the geometry of each bond type!!	
		if(INPUT.bdf_movie==1 and n_comm_neighbour==5 and nbonds==5)
		{
			static int count=0;
			if(count<1)
			{
				selected_atom[iat]=true;
				selected_atom[iat2]=true;
				for(int ii=0; ii < nadj_in[iat2]; ++ii)
				{
					const int iat3 = adj_index_in[iat2][ii];
					for(int jj=0; jj< nadj_in[iat]; ++jj)
					{
						const int iat4 = adj_index_in[iat][jj];
						if(iat3==iat4)
						{
							selected_atom[iat3]=true;
							//cout << iat3 << endl;
						}
					}
				}
			}
			++count; 
		}


		stringstream name;
		name << "1" << n_comm_neighbour << nbonds;


		for(int i=0; i<max_name; ++i)
		{
			if(HA_name[i]==name.str())
			{ 
				HA_percent[i]+=1.0; 
				break;
			}	
			else if(HA_name[i]=="none")
			{
				HA_name[i]=name.str();
				HA_percent[i]+=1.0;
				break;
			}
		}

//		cout << setw(10) << iat << " HA 1" << n_comm_neighbour << nbonds << endl;
//		ofs_running << setw(10) << iat << " HA 1" << n_comm_neighbour << nbonds << endl;
		

	}

	delete[] nnn;

	return;
}
