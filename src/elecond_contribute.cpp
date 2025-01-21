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
	double v11_all[NMAX],v12_all[NMAX],v22_all[NMAX], vweight[NMAX];
	for(int i=0;i<NMAX;i++)
	{
		L22_all[i]=0;
		L12_all[i]=0;
		sigma_all[i]=0;
		v11_all[i]=0;
		v12_all[i]=0;
		v22_all[i]=0;
		vweight[i]=0;
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
			// More electrons, smaller vmatrix
			// Since vmatrix means the transition probability，
			// more electrons, more choices, smaller probability for one choice.
			// Here we mutiple vmatrix by nele to normalize it.
			double vfact = pow(WF.factor,2)*INPUT.nele;
			//loop of band
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_sum, sigma_all[:NMAX], L12_all[:NMAX], L22_all[:NMAX], v11_all[:NMAX], v12_all[:NMAX], v22_all[:NMAX], vweight[:NMAX]) schedule(dynamic)
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
					double L11fact, L12fact, L22fact;
					double v11fact, v12fact, v22fact;
					double Eij_F=(energyj+energyi)/2-INPUT.fermiE;

					L11fact=factor2*docc*corr2;
					v11fact=vfact*corr2;

					L12fact=Eij_F*L11fact;
					L22fact=pow(Eij_F,2)*L11fact;
					v12fact=Eij_F*v11fact;
					v22fact=pow(Eij_F,2)*v11fact;
					if(INPUT.smear==0)
                	{
                	    std::cout<<"We do not support smear = 0 yet."<<std::endl;
						exit(1);
                	}
					else if(INPUT.smear==1||INPUT.smear==2)
					{
						double smearpre, invsmearpre;
					  	double pre;
						double st2;
						double eps2;
						if(INPUT.smear==1)
						{
							st2=-0.5/pow(st,2);
							pre=1/sqrt(2*pi)/st;
							double expminus = exp(pow(targew-w,2)*st2);
							double expplus = exp(pow(w+targew,2)*st2);
							smearpre = pre * (expminus - expplus);
							invsmearpre = pre * (expminus + expplus);
						}
						else
						{
							eps2=pow(eps,2);
							pre=1/pi*eps;
							double lminus = 1.0/(eps2+pow(targew-w,2));
							double lplus = 1.0/(eps2+pow(w+targew,2));
							smearpre = pre *(lminus - lplus);
							invsmearpre = pre *(lminus + lplus);
						}
						if(invsmearpre<1e-8) continue;
						int indexi = ((energyi - INPUT.fermiE) - INPUT.E_min)/INPUT.dE;
						int indexj = ((energyj - INPUT.fermiE) - INPUT.E_min)/INPUT.dE;
						double Lpre, vpre;
						if(INPUT.smearinvw)
						{
							Lpre = invsmearpre/w;
							vpre = smearpre/w*WF.wk;
						}
						else
						{
							Lpre = smearpre/targew;
							vpre = invsmearpre/targew*WF.wk;
						}
						if(indexi < NMAX && indexi >= 0)
						{
							sigma_all[indexi]+=Lpre*L11fact/2;
							L12_all[indexi]-=Lpre*L12fact/2;
							L22_all[indexi]+=Lpre*L22fact/2;
							v11_all[indexi]+=vpre*v11fact/2;
							v12_all[indexi]-=vpre*v12fact/2;
							v22_all[indexi]+=vpre*v22fact/2;
							vweight[indexi]+=vpre;
						}
						if(indexj < NMAX && indexj >= 0)
						{
							sigma_all[indexj]+=Lpre*L11fact/2;
							L12_all[indexj]-=Lpre*L12fact/2;
							L22_all[indexj]+=Lpre*L22fact/2;
							v11_all[indexj]+=vpre*v11fact/2;
							v12_all[indexj]-=vpre*v12fact/2;
							v22_all[indexj]+=vpre*v22fact/2;
							vweight[indexj]+=vpre;
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
	MPI_Allreduce(MPI_IN_PLACE, v11_all,NMAX,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, v12_all,NMAX,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, v22_all,NMAX,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, vweight,NMAX,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
#endif

	if(RANK == 0)
	{
		ofstream ofs1("contribute.txt");
		ofstream ofs2("vmatrix.txt");
		ofs1<<setw(10)<<"#E-mu (eV)"<<setw(20)<<"L11"<<setw(20)<<"L12/e"<<setw(20)<<"L22/e^2"<<endl;
		ofs2<<setw(10)<<"#E-mu (eV)"<<setw(20)<<"v11"<<setw(20)<<"v12/e"<<setw(20)<<"v22/e^2"<<setw(20)<<"weight"<<endl;
		double sum_sigma=0, sum_L12=0, sum_L22=0;
		for(int ie = 0 ; ie < NMAX; ++ie)
		{
			//contributions of Onsager coefficients
			ofs1<<setw(10)<<INPUT.E_min + (ie + 0.5) * INPUT.dE<<setw(20)<<sigma_all[ie] / nfolder/ INPUT.dE
			<<setw(20)<<L12_all[ie] / nfolder/ INPUT.dE<<setw(20)<<L22_all[ie] / nfolder/ INPUT.dE<<endl;
			sum_sigma += sigma_all[ie]/nfolder;
			sum_L12 += L12_all[ie]/nfolder;
			sum_L22 += L22_all[ie]/nfolder;

			//contributions of velocity matrix
			if(vweight[ie] > 0)
				ofs2<<setw(10)<<INPUT.E_min + (ie + 0.5) * INPUT.dE<<setw(20)<<v11_all[ie] / vweight[ie]/ INPUT.dE
			<<setw(20)<<v12_all[ie] / vweight[ie] / INPUT.dE<<setw(20)<<v22_all[ie] / vweight[ie] / INPUT.dE<<setw(20)<<vweight[ie] / nfolder<<endl;
			else
				ofs2<<setw(10)<<INPUT.E_min + (ie + 0.5) * INPUT.dE<<setw(20)<<0<<setw(20)<<0<<setw(20)<<0<<setw(20)<<0<<endl;
		}
		std::cout<<"sigma("<<targew<<") = "<<sum_sigma<<std::endl;
		std::cout<<"L12("<<targew<<") = "<<sum_L12<<std::endl;
		std::cout<<"L22("<<targew<<") = "<<sum_L22<<std::endl;
	}

    free(wdbuf);
	return;
}

	
