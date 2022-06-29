#include "cellFile.h"
#include "input.h"
#include "incremental_pdf.h"
#include "math.h"
#include "HBs.h"
#include "gfun.h"
#include "pdf_added.h"

incrementalPDF::incrementalPDF(){}
incrementalPDF::~incrementalPDF(){}
void incrementalPDF::Routine()
{
	TITLE("incrementalPDF", "Routine");
	cout << "Calculate the radial distribution functions g(r)." << endl;

	// (1) set 'delta r' in real space from input file.
	this->dr = INPUT.dr;
	assert(dr>0.0);

	// radius cutoff in real space, usually a/2,
	// where a is the lattice constant of a cubic box.
	this->rcut = INPUT.rcut;
	assert(rcut>0.0);

	// number of radial mesh grids.
	this->nmesh = int(rcut / dr) +  1;
	int nshell = INPUT.nshell;
	this->ntheta = int(INPUT.theta/INPUT.dtheta) + 1;
	this->nq = int(1/INPUT.dq)+1; // tetrahedral parameter q
	this-> tetraq_distr = new double[nq];
	for (int iq = 0; iq < nq; iq++)
	{
		tetraq_distr[iq] = 0.0;
	}
	ofs_running << "Number of shells considered: " << nshell << endl;
	// radial distribution function.
	double** gr;
	double** gr_tmp;
	gr = new double*[nshell];
	gr_tmp = new double*[nshell];
	this->theta_r_oo = new double**[nshell];
	this->tetraq_r_oo = new double**[nshell];
	int* total_nwater_in_shell = new int[nshell];
	double*** theta_r_oo_tmp = new double**[nshell];
	
	for(int i=0; i<nshell; ++i) 
	{
		gr[i] = new double[nmesh];
		gr_tmp[i] = new double[nmesh];
		total_nwater_in_shell[i] = 0;
		for(int j=0; j<nmesh; j++)
		{
			gr[i][j] = 0;
			gr_tmp[i][j] = 0;		
		}

		theta_r_oo[i] = new double*[ntheta];
		theta_r_oo_tmp[i] = new double*[ntheta];
		tetraq_r_oo[i] = new double*[nq];
		for (int itheta=0; itheta<ntheta; itheta++)
		{
			theta_r_oo[i][itheta] = new double[nmesh];
			theta_r_oo_tmp[i][itheta] = new double[nmesh];
			for(int imesh=0; imesh<nmesh; imesh++)
			{
				theta_r_oo[i][itheta][imesh] = 0.0;
				theta_r_oo_tmp[i][itheta][imesh] = 0.0;
			}
		}
		for(int iq = 0; iq<nq; iq++)
		{
			tetraq_r_oo[i][iq] = new double[nmesh];
			for(int imesh=0; imesh<nmesh; imesh++)
			{
				tetraq_r_oo[i][iq][imesh] = 0.0;
			}
		}
	}

	

	ofs_running << " dr = " << dr << " Angstrom" << endl;
	ofs_running << " rcut = " << rcut << " Angstrom" << endl;
	ofs_running << " nmesh = " << nmesh << endl;
	int count_geometry_number = 0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo) 
	{
		cout << "igeo=" << igeo << endl;
		CellFile cel;

		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;
		cout << "Succeeded" << endl;
		stringstream ss; ss << igeo;
		cel.file_name = ss.str();
		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) 
		{
			cel.clean();//qianrui add
			continue;
		}
		
		this->rho_ion = INPUT.natom/cel.volume;
		cout << "rho ion = " << rho_ion << endl;

		for(int ishell=0; ishell<nshell; ++ishell)
		{
			for(int imesh=0; imesh<nmesh; imesh++)
			{
				gr_tmp[ishell][imesh] = 0;
				for(int itheta=0; itheta<ntheta; itheta++)
				{
					theta_r_oo_tmp[ishell][itheta][imesh] = 0;
				}	
			}
		}

		//cout << "gr_tmp is set to 0." << endl;
		int ito = -1;
		int ith = -1;
		int itc = -1;
		int itcl = -1;
		int itna = -1;
		int it_ele1 = -1;
		int it_ele2 = -1;
		for(int it=0;it <INPUT.ntype; ++it)
		{
			if(cel.atom[it].id=="O") ito=it;
			else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
			else if(cel.atom[it].id=="C") itc=it;
			else if(cel.atom[it].id=="Cl") itcl=it;
			else if(cel.atom[it].id=="Na") itna=it;
			//cout << cel.atom[it].id << endl;
			if(INPUT.ele1 == cel.atom[it].id) it_ele1 = it;
			if(INPUT.ele2 == cel.atom[it].id) it_ele2 = it;
		}
		//cout << "ito = " << ito << "ith  = " << ith << endl;
	/*	for(int ia1=0; ia1<cel.atom[ito].na; ++ia1)
		{
			cout << cel.atom[ito].pos[ia1].x << " " << cel.atom[ito].pos[ia1].y << " " << cel.atom[ito].pos[ia1].z << endl;
		}
		for(int ia1=0; ia1<cel.atom[ith].na; ++ia1)
		{
			cout << cel.atom[ith].pos[ia1].x << " " << cel.atom[ith].pos[ia1].y << " " << cel.atom[ith].pos[ia1].z << endl;
		}*/
		Water *water;
		//cout << "water memory opened." << endl;
		// setup water molecules if necessary
			if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
		{
			water = new Water[cel.atom[ito].na];
			Water::nions = 0;
			//cout << "water memory = " << cel.atom[ito].na << endl;
			HBs::setup_water(cel, water);
			//cout << "Water setup done." << endl;
		}
		int io_of_ion = -1;
		if(INPUT.func == 1)
		// sort by distance
		{
			if(INPUT.system=="water")
			{
				this->distance_matrix = new double*[cel.atom[it_ele1].na];
				for (int ia1 = 0; ia1<cel.atom[it_ele1].na; ia1++)
				{
					this->distance_matrix[ia1] = new double[cel.atom[it_ele2].na];
					for (int ia2 = 0; ia2 < cel.atom[it_ele2].na; ia2++)
					{
						double dist = distance(cel.atom[it_ele1].pos[ia1], cel.atom[it_ele2].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						this->distance_matrix[ia1][ia2] = dist;	
					}
					bubble_sort(this->distance_matrix[ia1], cel.atom[it_ele2].na);
					// bubble sort each row of distance matrix
				} // end initializing distance matrix

				for (int ishell=0; ishell<nshell; ishell++)
				{
					int ia2 = -1;
					if(it_ele1==it_ele2) ia2 = ishell+1;
					// we don't wanna get the same atom involved in counting
					else ia2 = ishell;
					for(int ia1 = 0; ia1<cel.atom[it_ele1].na; ia1++)
					{
						int which = int(distance_matrix[ia1][ia2]/dr);
						gr_tmp[ishell][which]+=1.0;
					}
				}
				const double prec = 4.0/3.0*PI;
				for(int i=0; i<nmesh; ++i)
				{
					// volume
					double vv = prec*(pow((i+1)*dr,3)-pow(i*dr,3));
					const int n1 = cel.atom[it_ele1].na;
					const int n2 = cel.atom[it_ele2].na;
					for(int ishell=0; ishell<nshell; ishell++)
					{
						gr_tmp[ishell][i] = gr_tmp[ishell][i]*INPUT.natom/n1/n2/rho_ion/vv;
					}
				}

			}// end if water
			if(INPUT.system=="hydroxide" or INPUT.system=="hydronium")
			{	
				//cout << "start" << endl;
				if(INPUT.system=="hydronium" and INPUT.ele1=="H")
				{
					this->distance_matrix = new double*[7]();
					for (int ia1 = 0; ia1<7; ia1++)
					{
						this->distance_matrix[ia1] = new double[cel.atom[it_ele2].na]();
					}
				}
				else
				{
					this->distance_matrix = new double*[3]();
					this->distance_matrix[0] = new double[cel.atom[it_ele2].na]();	
					this->distance_matrix[1] = new double[cel.atom[it_ele2].na]();
					this->distance_matrix[2] = new double[cel.atom[it_ele2].na]();		
				}

	//			cout << "distance initialization done!" << endl;
				int na_ele1 = 0;
				for(int ia1=0; ia1<cel.atom[it_ele1].na; ++ia1)
				{

					// criterion 1: oxygen is in ion
					//cout << "na_ele1 = " << na_ele1 << endl;
					if(!PDF_ADDED::atom_in_ion(cel,it_ele1,ia1,io_of_ion)) continue;
					assert(io_of_ion>=0);
					//cout << "ion found!" << endl;
					for(int ia2 = 0; ia2<cel.atom[it_ele2].na; ia2++)
					{
						double dist = distance(cel.atom[it_ele1].pos[ia1], cel.atom[it_ele2].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						//cout << "dist = " << dist << endl;
						this->distance_matrix[na_ele1][ia2] = dist;
						
					}
					//cout << "finish diatance calc!" << endl;
					bubble_sort(this->distance_matrix[na_ele1], cel.atom[it_ele2].na);
					na_ele1++;
					//cout << "na_ele1 = " << na_ele1 << endl;
					cout << "finish bubble sorting!" << endl;
				}
				for (int ishell=0; ishell<nshell; ishell++)
				{
					int ia2 = -1;
					if(it_ele1==it_ele2) ia2 = ishell+1;
					else ia2 = ishell;
					for(int ia1 = 0; ia1<na_ele1; ia1++)
					{
						int which = int(this->distance_matrix[ia1][ia2]/dr);
						gr_tmp[ishell][which]+=1.0;
						//cout << distance_matrix[ia1][ia2] << endl;
						//cout << which << endl;
					}
				}
				//cout << "gr_tmp initialized" << endl;
				const double prec = 4.0/3.0*PI;
				for(int i=0; i<nmesh; ++i)
				{
					// volume
					double vv = prec*(pow((i+1)*dr,3)-pow(i*dr,3));
					const int n1 = 1;
					const int n2 = cel.atom[it_ele2].na;
					for(int ishell=0; ishell<nshell; ishell++)
					{
						gr_tmp[ishell][i] = gr_tmp[ishell][i]*INPUT.natom/n1/n2/rho_ion/vv;
					}
				}
			//cout << "gr " << gr_tmp[0][0] << endl;
			} // end if hydroxide or hydronium
		} // end if func == 1
		else if(INPUT.func == 2)
		{
			if(INPUT.system=="water")
			{
				//cout << "1" << endl;
				this->distance_matrix = new double*[cel.atom[it_ele1].na];
				this->index_matrix = new int*[cel.atom[it_ele1].na];
				for (int ia1 = 0; ia1<cel.atom[it_ele1].na; ia1++)
				{
					this->distance_matrix[ia1] = new double[cel.atom[ito].na];
					this->index_matrix[ia1] = new int[cel.atom[ito].na];
					for (int ia2 = 0; ia2 < cel.atom[ito].na; ia2++)
					{
						double dist = distance(cel.atom[it_ele1].pos[ia1], cel.atom[ito].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						this->distance_matrix[ia1][ia2] = dist;
						index_matrix[ia1][ia2] = ia2;	
					}
					//cout << "end initialization" << endl;
					bubble_sort2(this->distance_matrix[ia1], this->index_matrix[ia1], cel.atom[ito].na);
					// sort the order of indexed while sorting order of distances
				} // end initializing distance matrix
				//cout << "bubble sort done" << endl;
				for (int ishell=0; ishell<nshell; ishell++)
				{
					int ia2 = -1;
					if(it_ele1==it_ele2) ia2 = ishell+1;
					else ia2 = ishell;
					for(int ia1 = 0; ia1<cel.atom[it_ele1].na; ia1++)
					{
						for(int iah = 0; iah<water[index_matrix[ia1][ia2]].nH; iah++)
						{
							int iih = water[index_matrix[ia1][ia2]].indexH[iah];
							double dist1 = distance(cel.atom[it_ele1].pos[ia1], cel.atom[ith].pos[iih], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
							int which = int(dist1/dr);
							gr_tmp[ishell][which]+=1.0;
						}//end for iah	
					}//end for ia1
				}//end for ishell
				const double prec = 4.0/3.0*PI;
				for(int i=0; i<nmesh; ++i)
				{
					// volume
					double vv = prec*(pow((i+1)*dr,3)-pow(i*dr,3));
					const int n1 = cel.atom[it_ele1].na;
					const int n2 = cel.atom[it_ele2].na;
					for(int ishell=0; ishell<nshell; ishell++)
					{
						gr_tmp[ishell][i] = gr_tmp[ishell][i]*INPUT.natom/n1/n2/rho_ion/vv;
					}
				}
				for(int ia1 = 0; ia1<cel.atom[it_ele1].na; ia1++)
				{
					delete[] this->distance_matrix[ia1];
					delete[] this->index_matrix[ia1];
				}
				delete[] this->distance_matrix;
				delete[] this->index_matrix;
			}//end if water
			else if(INPUT.system=="hydroxide" or INPUT.system=="hydronium")
			{	
				cout << "start" << endl;
				this->distance_matrix = new double*[1];
				this->index_matrix = new int*[1];
				distance_matrix[0] = new double[cel.atom[ito].na];
				index_matrix[0] = new int[cel.atom[ito].na];
	//			cout << "distance initialization done!" << endl;
				int na_ele1 = 0;
				if(Water::nions != 1) continue;
				int stop = -1;
				int ia_of_ion = -1;
				for(int ia1=0; ia1<cel.atom[ito].na; ++ia1)
				{
					if((water[ia1].nH==3 and INPUT.system=="hydronium") or (water[ia1].nH==1 and INPUT.system=="hydroxide"))
					{
						if((INPUT.nacc==3 and water[ia1].naccept!=3) or (INPUT.nacc==4 and water[ia1].naccept!=4) or (INPUT.nacc==5 and water[ia1].naccept!=5))
						{
							stop = 1;
							continue;
						}
						for (int ia2 = 0; ia2 < cel.atom[ito].na; ia2++)
						{
							double dist = distance(cel.atom[ito].pos[ia1], cel.atom[ito].pos[ia2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
							this->distance_matrix[0][ia2] = dist;
							index_matrix[0][ia2] = ia2;	
						}
						//cout << "end initialization" << endl;
						bubble_sort2(this->distance_matrix[0], this->index_matrix[0], cel.atom[ito].na);
						ia_of_ion = ia1;
						break; // no need to loop any more
					} //end if ion
				} //end ia1
				if(stop==1) continue;// if the input condition is not satisfied, withdraw from the loop
				//cout << water[ia_of_ion].naccept << endl;
				for (int ia2=1; ia2<nshell; ia2++)
				{
					for(int iah = 0; iah<water[index_matrix[0][ia2]].nH; iah++)
						{
							int iih = water[index_matrix[0][ia2]].indexH[iah];
							double dist1 = distance(cel.atom[ito].pos[ia_of_ion], cel.atom[ith].pos[iih], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
							int which = int(dist1/dr);
							gr_tmp[ia2-1][which]+=1.0;
						}//end for iah	
				}//end for ishell
				cout << "gr_tmp initialized" << endl;
				const double prec = 4.0/3.0*PI;
				for(int i=0; i<nmesh; ++i)
				{
					// volume
					double vv = prec*(pow((i+1)*dr,3)-pow(i*dr,3));
					const int n1 = 1;
					const int n2 = cel.atom[it_ele2].na;
					for(int ishell=0; ishell<nshell; ishell++)
					{
						gr_tmp[ishell][i] = gr_tmp[ishell][i]*INPUT.natom/n1/n2/rho_ion/vv;
					}
				}
				delete[] this->distance_matrix[0];
				delete[] this->index_matrix[0];

				delete[] this-> distance_matrix;
				delete[] this-> index_matrix;
			} // end if hydroxide or hydronium
		} // end func=2
		else if(INPUT.func == 3)
		// topological neighbor sort (order of topological neighbor is described in 15PCCP-Imoto)
		{
			this->topological_neighbor_PDF(cel, water, gr_tmp, total_nwater_in_shell, nshell, ito, ith);
		//	cout << "func = 3 calculation done." << endl;
		}// end func==3
		for(int i = 0; i<nmesh; i++)
		{
			for(int ishell = 0; ishell<nshell; ishell++)
			{
				gr[ishell][i]+=gr_tmp[ishell][i];
				/*for (int itheta=0; itheta<ntheta; itheta++)
				{
					theta_r_oo[ishell][itheta][i] += theta_r_oo_tmp[ishell][itheta][i];
				}*/
			}
		}
		

		count_geometry_number++;
		cel.clean();
		
		cout << "igeo " << igeo << endl;
		ofs_running << "count_geometry_number=" << count_geometry_number << endl;
		if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide")
		{
			delete[] water;
		}
	}// igeo
	assert(count_geometry_number>0);
    ofs_running << " count_geometry_number = " << count_geometry_number << endl;

    	if(count_geometry_number>0)
    	{
    		for(int ir=0; ir<nmesh; ++ir)
    	    {
				for(int ishell=0; ishell<nshell; ishell++)
				{
					gr[ishell][ir] /= count_geometry_number;
					for (int itheta=0; itheta<ntheta; itheta++)
					{
						theta_r_oo[ishell][itheta][ir] /= count_geometry_number;
					}
					if (INPUT.func_b == 2)
					{
						for (int iq = 0; iq < nq; iq++)
						{
							tetraq_r_oo[ishell][iq][ir] /= count_geometry_number;
						}
					}
				}
    	    }
    	}
   	for (int ishell = 0; ishell<nshell; ishell++)
	{
		
		double sum = 0;
		for (int itheta=0; itheta<ntheta; itheta++)
		{
			for(int imesh=0; imesh<nmesh-1; imesh++)
			{
				sum += INPUT.dtheta*INPUT.dr*theta_r_oo[ishell][itheta][imesh];
			}
		}
		
		ofstream ofs_theta_r_oo("theta_r_oo_shell" + to_string(ishell)+".txt");
		for (int ir=0; ir<nmesh-1; ir++)
		{
			ofs_theta_r_oo << "  " << ir*INPUT.dr+0.5*INPUT.dr;
		}
		ofs_theta_r_oo << endl;
		for (int itheta=0; itheta<ntheta; itheta++)
		{
			ofs_theta_r_oo << itheta*INPUT.dtheta+0.5*INPUT.dtheta << " ";
			for(int imesh=0; imesh<nmesh-1; imesh++)
			{
				double rho = INPUT.natom1/INPUT.celldm1/INPUT.celldm2/INPUT.celldm3;
            	double norm_factor = rho*double(count_geometry_number)*INPUT.dtheta*INPUT.dr*2*PI*(imesh+0.5)*(imesh+0.5)*INPUT.dr*INPUT.dr*sin((itheta+0.5)*INPUT.dtheta/180*PI);
				ofs_theta_r_oo << theta_r_oo[ishell][itheta][imesh]/norm_factor << " " ;
			}
			ofs_theta_r_oo << endl;
		}
		ofs_theta_r_oo.close();

		if (INPUT.func_b == 2)
		{
			double sumq = 0;
			for (int iq=0; iq<nq; iq++)
			{
				for(int imesh=0; imesh<nmesh-1; imesh++)
				{
					sum += INPUT.dq*INPUT.dr*tetraq_r_oo[ishell][iq][imesh];
				}
			}
			ofstream ofs_tetraq_r_oo("tetraq_r_oo_shell" + to_string(ishell)+".txt");
			for (int ir=0; ir<nmesh-1; ir++)
			{
				ofs_tetraq_r_oo << "  " << ir*INPUT.dr+0.5*INPUT.dr;
			}
			ofs_tetraq_r_oo << endl;
			for (int iq=0; iq<nq; iq++)
			{
				ofs_tetraq_r_oo << iq*INPUT.dq+0.5*INPUT.dq << " ";
				for(int imesh=0; imesh<nmesh-1; imesh++)
				{
					ofs_tetraq_r_oo << tetraq_r_oo[ishell][iq][imesh]/sumq << " " ;
				}
				ofs_tetraq_r_oo << endl;
			}
			ofs_tetraq_r_oo.close();
		}
   }
	ofstream ofs(INPUT.geo_out.c_str());
	for(int imesh = 0; imesh<nmesh-1; imesh++)
	{	
		ofs << (imesh+0.5)*dr << " ";
		for(int ishell = 0; ishell < nshell; ishell++)
		{
			ofs << gr[ishell][imesh] << " ";	
		}
		ofs << endl;
	}
	ofs.close();
	
	ofs_running << "shell    nwater" << endl;
	for(int ishell=0; ishell<nshell; ++ishell)
	{ 
		ofs_running << ishell << "    " << total_nwater_in_shell[ishell]/count_geometry_number << endl;
	}

	if (INPUT.func_b == 2)
	{
		ofstream ofs_tetraq("tetraq_distribution.txt");
		double sum = 0;
		for (int iq = 0; iq < nq; iq++)
		{
			tetraq_distr[iq] /= count_geometry_number;
			sum += tetraq_distr[iq]*INPUT.dq;
		}
		double average = 0;
		for (int iq = 0; iq < nq; iq++)
		{
			ofs_tetraq << iq*INPUT.dq + 0.5*INPUT.dq << " " << tetraq_distr[iq] / sum << endl;
			average += (iq*INPUT.dq + 0.5*INPUT.dq)*tetraq_distr[iq]/sum*INPUT.dq;
		}
		ofs_running << "average tetrahedral q = " << average << endl;
		ofs_tetraq.close();
	}

	//	return;
}

void incrementalPDF::bubble_sort(double* list, const int &length)
{
	int i, j;
	for (i = 0; i<length-1; i++)
	{
		for (j = 0; j<length-1-i; j++)
		{
			if(list[j]>list[j+1])
			{
				double tmp = list[j];
				list[j] = list[j+1];
				list[j+1] = tmp;	
			}
		}
	}
}

void incrementalPDF::bubble_sort2(double* list, int* index_list, const int &length)
{
	int i, j;
	for (i = 0; i<length-1; i++)
	{
		for (j = 0; j<length-1-i; j++)
		{
			if(list[j]>list[j+1])
			{
				double tmp = list[j];
				list[j] = list[j+1];
				list[j+1] = tmp;

				int index_tmp = index_list[j];
				index_list[j] = index_list[j+1];
				index_list[j+1] = index_tmp;	
			}
		}
	}
}

void incrementalPDF::topological_neighbor_PDF(CellFile &cel, Water* &water, double** &gr_tmp, int* &total_nwater_in_shell, int nshell, int ito, int ith)
{
	//cout << "func = " << INPUT.func << endl;

	if(INPUT.system == "water")
	{
		// Step 1: Initialize essential int lists.
		int* included_water = new int[cel.atom[ito].na];
		int nincluded_water = 0;

		int* nwater_in_shell = new int[nshell];
		int** water_in_shell = new int*[nshell];
		for(int ishell=0; ishell<nshell; ++ishell)
		{
			nwater_in_shell[ishell] = 0;
			water_in_shell[ishell] = new int[cel.atom[ito].na];
		}
		for(int ia1=0; ia1<cel.atom[ito].na; ++ia1)
		{
			double tetraq = 1;
			if (INPUT.func_b == 2)
			{
				// calculate q parameter
				// q parameter is for every central water molecule, so can be calculated before shell calculation
				double* distance_arr = new double [cel.atom[ito].na];
				int* index_arr = new int [cel.atom[ito].na];
				for (int ia = 0; ia <cel.atom[ito].na; ++ia)
				{
					distance_arr[ia] = 0;
					index_arr[ia] = ia;
					if(ia != ia1)
					{
						distance_arr[ia] = distance(cel.atom[ito].pos[ia], cel.atom[ito].pos[ia1], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
					}
				}
				bubble_sort2(distance_arr, index_arr, cel.atom[ito].na);
				assert(index_arr[0] == ia1);
				for (int iaa=1; iaa<4; iaa++)
				{
					for (int iaa2=iaa+1; iaa2 < 5; iaa2++)
					{
						int index1 = index_arr[iaa];
						int index2 = index_arr[iaa2];
						double dist_c = distance(cel.atom[ito].pos[index1], cel.atom[ito].pos[index2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						double cos = (distance_arr[iaa]*distance_arr[iaa] + distance_arr[iaa2]*distance_arr[iaa2] - dist_c*dist_c)/2/distance_arr[iaa]/distance_arr[iaa2];
						//cout << distance_arr[iaa] << " " << distance_arr[iaa2] << " " << dist_c << endl;
						//cout << "iaa = " << iaa << ", iaa2 = " << iaa2 << ", cos = " << cos << endl;
						tetraq -= 3*(cos+1/3)*(cos+1/3)/8;
					}
				}
				//cout << "ia = " << ia1 << ", tetraq = " << tetraq << endl;
				if (tetraq < 0)
				{
					if (tetraq > -0.1)
					{
						tetraq = 0;
					}
					else
					{
						continue;
						//cout << "tetraq = " << tetraq << endl;
						//exit(0);
					}
				}
				if (tetraq > 1)
				{
					if (tetraq < 1.1)
					{
						tetraq = 1;
					}
					else
					{
						cout << "tetraq = " << tetraq << endl;
						exit(0);
					}
				}
				int which_q = (int) (tetraq/INPUT.dq);
				tetraq_distr[which_q]++;
				delete[] distance_arr;
				delete[] index_arr;
			}// if func_b == 2 calculate tetraq
			nincluded_water = 0;
			for(int ia=0; ia<cel.atom[ito].na; ++ia)
			{
				included_water[ia] = -1;
			}
			for(int ishell=0; ishell<nshell; ++ishell)
			{
				for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
				{
					water_in_shell[ishell][ia2] = -1;
				}// ia2
			}//ishell
			included_water[nincluded_water] = ia1;
			++nincluded_water;
			int nwater_in_last_shell;
			for(int ishell=0; ishell<nshell; ++ishell)
			{
				//cout << "new iter ishell = " << ishell << " begins." << endl;
				if(ishell==0)
				{
					//cout << "0 shell" << endl;
					nwater_in_last_shell = 0;
					for(int ia2=0; ia2<water[ia1].naccept;++ia2)
					{
						water_in_shell[ishell][ia2] = water[ia1].acceptO[ia2];
						
						included_water[nincluded_water] = water[ia1].acceptO[ia2];
						//cout << ia2 << " " << water_in_shell[ishell][ia2] << endl;
						++nincluded_water;
					}
					for(int ia2=0; ia2<water[ia1].ndonate; ++ia2)
					{
						int water_index = water[ia1].donateO[ia2];
						bool included = false;
						for(int ia_check=0; ia_check<water[ia1].naccept; ++ia_check)
						{
							if( water[ia1].acceptO[ia_check] == water_index )
							{
								included = true;
								continue;
							}
						}
						if(!included)
						{
							//water_in_shell[ishell][ia2+water[ia1].naccept] = water_index;
							water_in_shell[ishell][nincluded_water-1] = water_index;
							//cout << nincluded_water-1 << " " << water_in_shell[ishell][nincluded_water-1] << endl;
							included_water[nincluded_water] = water_index;
							++nincluded_water;
						}
					}
					//cout << ia1 << " naccept = " << water[ia1].naccept << " ndonate = " << water[ia1].ndonate << endl;
					nwater_in_last_shell = nincluded_water-1;
					nwater_in_shell[ishell] += nwater_in_last_shell;
					//cout << "water " << ia1 << " has " << nwater_in_shell[ishell] << " water in shell " << ishell << endl;
					/*
					cout << "nwater_in_last_shell = " << nwater_in_last_shell << endl;
					cout << "accept ";
					for(int ia_check=0; ia_check<water[ia1].naccept; ++ia_check)
					{
						cout << water[ia1].acceptO[ia_check] << " " ;
					}
					cout << endl << "donate ";
					for(int ia_check=0; ia_check<water[ia1].ndonate; ++ia_check)
					{
						cout << water[ia1].donateO[ia_check] << " " ;
					}
					cout << endl;*/
				}
				else if(ishell < nshell-1)
				{	
					//cout << "ishell = " << ishell << endl;
					//cout << "nwater_in_last_shell="<<nwater_in_last_shell << endl;
					int nwater_in_this_shell = 0;
					//cout << "nwater_in_last_shell = " << nwater_in_last_shell << endl;
					for(int ia2=0; ia2<nwater_in_last_shell; ++ia2)
					{
						int water_index1 = water_in_shell[ishell-1][ia2];
						//cout << "water_index1 = " << water_index1 << endl;
						//cout << ia1 << " " << ishell << " " << ia2 << endl;
						/*
						cout << "water ia2 = " << ia2 << " water_index1 = " << water_index1 <<  " has " << water[water_index1].naccept << " accepted waters." << endl;
						cout << "water ia2 = " << ia2 << " water_index1 = " << water_index1 <<  " has " << water[water_index1].ndonate << " donated waters." << endl;
						*/
						for(int ia3=0; ia3<water[water_index1].naccept; ++ia3)
						{
							int water_index2 = water[water_index1].acceptO[ia3];
							bool included = false;
							for(int ia_check=0; ia_check<nincluded_water; ++ia_check)
							{
								if(included_water[ia_check] == water_index2)
								{
									included = true;
									continue;
								}
							}
							if(!included)
							{
								water_in_shell[ishell][nwater_in_this_shell] = water_index2;
								++nwater_in_this_shell;
								included_water[nincluded_water] = water_index2;
								++nincluded_water;
								if(nincluded_water == cel.atom[ito].na)
								{
									cout << "Waters within " << ishell << "th shell accept are enough to cover all waters." << endl;
								}
							}
						}// ia3 accept
						//cout << "ia1 = " << ia1 << ", ishell = " << ishell << " accept done." << endl;
						for(int ia3=0; ia3<water[water_index1].ndonate; ++ia3)
						{
							int water_index2 = water[water_index1].donateO[ia3];
							bool included = false;
							for(int ia_check=0; ia_check<nincluded_water; ++ia_check)
							{
								if(included_water[ia_check] == water_index2)
								{
									included = true;
								}
							}
							if(!included)
							{
								water_in_shell[ishell][nwater_in_this_shell] = water_index2;
								++nwater_in_this_shell;
								included_water[nincluded_water] = water_index2;
								++nincluded_water;
								if(nincluded_water == cel.atom[ito].na)
								{
									cout << "Waters within " << ishell << "th shell donate are enough to cover all waters." << endl;
									cout << "nwater_in_last_shell = " << nwater_in_last_shell << endl;
									cout << "water_in_this_shell = " << nwater_in_this_shell << endl;
									//exit(0);
								}
							}
						}// ia3 donate
					}// ia2
					nwater_in_shell[ishell] += nwater_in_this_shell;
					nwater_in_last_shell = nwater_in_this_shell;
					//cout << "water " << ia1 << " has " << nwater_in_this_shell << " water in shell " << ishell << endl;
					if(nincluded_water >= cel.atom[ito].na)
					{
						cout << "Waters within " << nshell << " shells are enough to cover all waters." << endl;
						//exit(0); 
					}
				}// if ishell < nshell
				else
				{
					//cout << ia1 << " come to final shell" << endl;
					assert(ishell == nshell-1);
					int nwater_out_of_shell = 0;
					for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
					{
						int out = true;
						for(int ia3=0; ia3<nincluded_water; ++ia3)
						{
							if(ia2 == included_water[ia3])
							{
								out = false;
							}
						}
						if(out)
						{
							water_in_shell[nshell-1][nwater_out_of_shell] = ia2;
							++nwater_out_of_shell;
						}
					}// ia2
					if(nwater_out_of_shell+nincluded_water!=cel.atom[ito].na)
					{
						cout << "nwater_out_of_shell = " << nwater_out_of_shell << endl;
						cout << "nincluded_water = " << nincluded_water << endl;
						cout << "ia1 = " << ia1 << endl;
						int sum = ia1;
						for(int ishell=0; ishell<nshell; ++ishell)
						{
							cout << ishell << endl;
							for(int ia_check=0; ia_check<nwater_in_shell[ishell]; ++ia_check)
							{
								cout << water_in_shell[ishell][ia_check] << " ";
								if(water_in_shell[ishell][ia_check]>=0 and water_in_shell[ishell][ia_check]<cel.atom[ito].na)
								{sum += water_in_shell[ishell][ia_check];}
							}
							cout << endl;
							cout << "sum = " << sum << endl;
						}
						exit(0);
					}
					nwater_in_shell[nshell-1] += nwater_out_of_shell;
				} // last shell
			}//ishell
			
			// calculate gr_tmp
			
			for(int ishell=0; ishell<nshell; ++ishell)
			{
				//cout << "ishell = " << ishell << endl;
				for(int ia2=0; ia2<cel.atom[ito].na; ++ia2)
				{
					int water_index = water_in_shell[ishell][ia2];
					if(water_index != -1 and water_index < cel.atom[ito].na)
					{
						//cout << water_index << " ";
						//cout << ishell << " " << water_index << endl;
						double dist1 = distance(cel.atom[ito].pos[ia1], cel.atom[ito].pos[water_index], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
						int which = int(dist1/dr);
						if(which < nmesh)
						{
							gr_tmp[ishell][which]+=1.0;
						}
						// renxi 20210726
						int Hindex0=water[ia1].indexH[0];
						int Hindex1=water[ia1].indexH[1];
						//cout << ia1 << " " << Hindex0 << " " << Hindex1 << endl;
						double angle0 = HBs::angle(cel, cel.atom[ith].pos[Hindex0], cel.atom[ito].pos[ia1], cel.atom[ito].pos[water_index]);
						double angle1 = HBs::angle(cel, cel.atom[ith].pos[Hindex1], cel.atom[ito].pos[ia1], cel.atom[ito].pos[water_index]);
						int which_angle0 = int(angle0/INPUT.dtheta);
						int which_angle1 = int(angle1/INPUT.dtheta);
						//cout << angle0 << " " << angle1 << endl;
						this->theta_r_oo[ishell][which_angle0][which]++;
						this->theta_r_oo[ishell][which_angle1][which]++;
						//cout << this->theta_r_oo[ishell][which_angle0][which] << endl;

						// renxi 20210730 malong yyds
						// tetraq
						if (INPUT.func_b == 2)
						{
							assert(tetraq >= 0 and tetraq <= 1);
							int which_q = (int) (tetraq/INPUT.dq);
							this->tetraq_r_oo[ishell][which_q][which]++;
						}
					}
					else
					{
						break;
					}
				}
				//cout << endl;
			}// ishell for output
			/*
			cout << ia1 << " gr_tmp = " << endl;
			for(int ishell = 0; ishell<nshell; ++ishell)
			{
				int sum = 0;
				for(int imesh=0; imesh<nmesh; ++imesh)
				{
					sum += gr_tmp[ishell][imesh];
				}
				cout << "ishell = " << ishell << " sum = " << sum << endl;
			}*/
			
		}//ia1
		//cout << "ia1 done." << endl;
		for(int ishell=0; ishell<nshell; ++ishell)
		{
			nwater_in_shell[ishell] /= cel.atom[ito].na;
			total_nwater_in_shell[ishell] += nwater_in_shell[ishell];
			//cout << "ishell = " << ishell << " waters_in_shell[ishell] = " << nwater_in_shell[ishell] << endl;
		
		//cout << "number of water calculation done." << endl;
			const double prec = 4.0/3.0*PI;
		
			for(int i=0; i<nmesh; ++i)
			{
				// volume
				double vv = prec*(pow((i+1)*dr,3)-pow(i*dr,3));
				const int n1 = cel.atom[ito].na;
					/*if(gr_tmp[ishell][i] != 0)
					{
						cout << "ishell = " << ishell << endl;
						cout << "i = " << i << endl;
						cout << "INPUT.natom = " << INPUT.natom << endl;
						cout << "n1 = " << n1 << endl;
						cout << "rho_ion = " << rho_ion << endl;
						cout << "vv = " << vv << endl; 
						cout << "gr_tmp[ishell][i] = " << gr_tmp[ishell][i] << endl;
						cout << "gr_tmp[ishell][i] = " << gr_tmp[ishell][i]*INPUT.natom/n1/n1/rho_ion/vv << endl;
					}*/
					gr_tmp[ishell][i] = gr_tmp[ishell][i]*INPUT.natom/n1/n1/rho_ion/vv;
				
			}// imesh
		}// ishell
		/*
		for(int ishell = 0; ishell<nshell; ++ishell)
		{
			int sum = 0;
			for(int imesh=0; imesh<nmesh; ++imesh)
			{
				sum += gr_tmp[ishell][imesh];
			}
			cout << "Final gr_tmp = " << endl;
			cout << "ishell = " << ishell << " sum = " << sum << endl;
		}
		
		for(int ishell=0; ishell<nshell; ++ishell)
		{
			cout << ishell << endl;
			for(int i=0; i<nmesh; ++i)
			{
				cout << gr_tmp[ishell][i] << " ";
			}
			cout << endl;
		}
		*/
		//cout << sizeof(gr_tmp)/sizeof(gr_tmp[0]) << " " << sizeof(gr_tmp[0])/sizeof(gr_tmp[0][0]) << endl;



		// Clear all memories
		delete[] included_water;
		delete[] nwater_in_shell;
		for(int ishell=0; ishell<nshell; ++ishell)
		{
			delete[] water_in_shell[ishell];
		}
		delete[] water_in_shell;
		//cout << "Function done." << endl;


	}
	return;
}
