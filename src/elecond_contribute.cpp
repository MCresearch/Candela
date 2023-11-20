#include "input.h"
#include "const.h"
#include "gfun.h"
#include "elecond_contribute.h"
#ifdef __MPI
#include "mpi.h"
#endif
#include <unistd.h>

void Elecond_contribute::Routine()
{
	this->cal_contribute();
	return;
}

void Elecond_contribute::cal_contribute()
{
	//get current work directory 
	char *wdbuf;
	wdbuf=getcwd(nullptr,0);
	//cout<<wdbuf<<endl;

	//Init
	const double pi=M_PI;
	double fwhm = INPUT.fwhm[0];
	double st = fwhm/TWOSQRT2LN2, eps = fwhm/2;
	double factor = 2*pi/3*pow(P_HBAR,3)/(pow(P_ME,2)*pow(P_BOHR*1e-10,5));
	int allcount=-1;

	int nfolder=INPUT.nscf;//number of folders containing different scf results.
	bool multi=true;
	double tot_sum=0;
	if(nfolder<=0)
	{ 
		nfolder=1;
		multi=false;
	}
	
	double targew = INPUT.target_w;
	
	int NMAX = (INPUT.E_max - INPUT.E_min)/INPUT.dE + 1;
	assert(NMAX > 0);
	double sigma_all[NMAX],L22_all[NMAX],L12_all[NMAX];
	for(int i=0;i<NMAX;i++)
	{
		L22_all[i]=0;
		L12_all[i]=0;
		sigma_all[i]=0;
	}

	if(multi)
	{
		string command="test -e "+INPUT.multi_directory;
		int exist=system(command.c_str());
		if(exist!=0)
		{
			cout<<"No "<<INPUT.multi_directory<<" exists."<<endl;
			exit(0);
		}
		chdir(INPUT.multi_directory.c_str());
	}
	
	//loop of snapshot(scf)
	for(int ifolder=0;ifolder<nfolder;ifolder++)
	{
		if(multi)
		{
			int ifoldername=ifolder+1;
			string ifolderstr=int2str(ifoldername);
			string command="test -e "+ifolderstr;
			int exist=system(command.c_str());
			if(exist!=0)
			{
				cout<<"No "<<ifolderstr<<" exists."<<endl;
				exit(0);
			}
			chdir(ifolderstr.c_str());
		}

		Wavefunc WF;
		WfRead wfr(WF);
		wfr.Init();
		double sum_factor=2*P_ME*INPUT.vol*pow(P_BOHR,3)/M_PI/P_QE/INPUT.nele/1e30/P_HBAR;
		int nk=INPUT.nkpoint;
		assert(nk>0);

		//kpioint loop
		for(int ik=0;ik<nk;ik++)
		{
			allcount++;
#ifdef __MPI
        	if((allcount)%NPROC!=RANK)
			{
				wfr.ignore(ik);
				wfr.clean();
				continue;
			}
#endif
			wfr.readWF(ik);
			int nband=WF.nband;
			cout<<"scf "<<ifolder<<" ; kpoint "<<ik+1<<endl;
			//WF.print(0);

			const int nbb = (nband-1) * nband / 2;
			double *vmatrix = new double [nbb];
			if(INPUT.readvmatrix)
			{
				wfr.readvmatrix(ik, vmatrix);
			}
			else
			{
#ifdef _OPENMP
#pragma omp parallel for
#endif	
				for(int ib=0;ib<nband;ib++)
				{
					WF.checknorm(ik,ib);//check if wavefunction is normalized to 1.
				}
				wfr.calvmatrix(vmatrix);
			}
			cout<<"nband: "<<nband<<endl;
			double factor2=factor*pow(WF.factor,2)/INPUT.vol;
			//loop of band
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_sum, sigma_all[:NMAX], L12_all[:NMAX], L22_all[:NMAX]) schedule(dynamic)
#endif	
			for(int ib=0;ib<nband;++ib)
			{
				for(int jb=ib+1;jb<nband;++jb)
	    		{
					int ijb = nband*ib - ib*(ib+1)/2 + (jb - ib - 1);
					double energyj=WF.eigE[jb];
					double energyi=WF.eigE[ib];
					double w=energyj-energyi;
					double docc=WF.occ[ib]-WF.occ[jb];
					if(docc<=0.0000000001 || w < 1e-5) continue;

					double corr2 = vmatrix[ijb];

//--------------------------------------------------------------------------------------------------
//Smearing Part
//Frequency-dependent functions are (or not (0)) broadened by Gaussian (0) or Lorentz (2) function.
//--------------------------------------------------------------------------------------------------
					double L11fact;
					double Eij_F=(energyj+energyi)/2-INPUT.fermiE;

					if(INPUT.smearinvw)
					{
						L11fact=factor2/w*docc*corr2;
					}
					else
					{
						L11fact=factor2*docc*corr2;
					}

					double L12fact=Eij_F*L11fact;
					double L22fact=pow(Eij_F,2)*L11fact;
					if(INPUT.smear==0)
                	{
                	    std::cout<<"We do not support smear = 0 yet."<<std::endl;
						exit(1);
                	}
					else if(INPUT.smear==1||INPUT.smear==2)
					{
					  	double smearf;
						double smearpre;
					  	double pre;
						double st2;
						double eps2;
						if(INPUT.smear==1)
						{
							st2=-0.5/pow(st,2);
							pre=1/sqrt(2*pi)/st;
							if(INPUT.smearinvw)
							{
								smearf=exp(pow(targew-w,2)*st2) + exp(pow(w+targew,2)*st2);
								smearpre = pre * smearf;
							}
							else
							{
								smearf=exp(pow(targew-w,2)*st2) - exp(pow(w+targew,2)*st2);
								smearpre = pre/targew * smearf;
							}
						}
						else
						{
							eps2=pow(eps,2);
							pre=1/pi*eps;
							if(INPUT.smearinvw)
							{
								smearf = 1/(eps2+pow(targew-w,2)) + 1/(eps2+pow(w+targew,2));
								smearpre = pre * smearf;
							}
							else
							{
								smearf = 1/(eps2+pow(targew-w,2)) - 1/(eps2+pow(w+targew,2));
								smearpre = pre/targew * smearf;
							}
						}
						if(smearf<1e-6) continue;
						int indexi = ((energyi - INPUT.fermiE) - INPUT.E_min)/INPUT.dE;
						int indexj = ((energyj - INPUT.fermiE) - INPUT.E_min)/INPUT.dE;
						if(indexi < NMAX && indexi >= 0)
						{
							sigma_all[indexi]+=smearpre*L11fact/2;
							L12_all[indexi]-=smearpre*L12fact/2;
							L22_all[indexi]+=smearpre*L22fact/2;
						}
						if(indexj < NMAX && indexj >= 0)
						{
							sigma_all[indexj]+=smearpre*L11fact/2;
							L12_all[indexj]-=smearpre*L12fact/2;
							L22_all[indexj]+=smearpre*L22fact/2;
						}
					}
				}//jb
	    	}//ib
			wfr.clean();
			delete[] vmatrix;
		}//ik
		wfr.cleanclass();
		if(multi)
		{
			chdir("../");
		}
	}//ifolder

	if(multi)
	{
		chdir(wdbuf);
	}
#ifdef __MPI
	MPI_Allreduce(MPI_IN_PLACE, sigma_all,NMAX,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, L12_all,NMAX,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, L22_all,NMAX,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

	if(RANK == 0)
	{
		ofstream ofs1("contribute.txt");
		ofs1<<setw(10)<<"#E-mu (eV)"<<setw(20)<<"S(meV)^-1"<<setw(20)<<"A(meV)^-1"<<setw(20)<<"W(mKeV)^-1"<<endl;
		double sum_sigma=0, sum_L12=0, sum_L22=0;
		for(int ie = 0 ; ie < NMAX; ++ie)
		{
			ofs1<<setw(10)<<INPUT.E_min + (ie + 0.5) * INPUT.dE<<setw(20)<<sigma_all[ie] / nfolder/ INPUT.dE
			<<setw(20)<<L12_all[ie] / nfolder/ INPUT.dE<<setw(20)<<L22_all[ie] / nfolder/ INPUT.dE<<endl;
			sum_sigma += sigma_all[ie]/nfolder;
			sum_L12 += L12_all[ie]/nfolder;
			sum_L22 += L22_all[ie]/nfolder;
		}
		std::cout<<"sigma("<<targew<<") = "<<sum_sigma<<std::endl;
		std::cout<<"L12("<<targew<<") = "<<sum_L12<<std::endl;
		std::cout<<"L22("<<targew<<") = "<<sum_L22<<std::endl;
	}

    free(wdbuf);
	return;
}

	
