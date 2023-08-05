#include "input.h"
#include "const.h"
#include "gfun.h"
#include "ele_conductivity.h"
#include <unistd.h>

//sigma is electric conductivity and kappa is thermal conductivity.
//L_{mn}(\omega)=(-1)^{m+n}\frac{2\pi e^2\hbar^2}{3m_e^2\omega\Omega}\\
//\times\sum_{i,j,\alpha,\mathbf{k}}W(\mathbf{k})\left(\frac{\epsilon_{i,\mathbf{k}}+\epsilon_{j,\mathbf{k}}}{2}-\mu\right)^{m+n-2}|\langle\Psi_{i,\mathbf{k}}|\nabla_\alpha|\Psi_{j,\mathbf{k}}\rangle|^2\\
//\times[f(\epsilon_{i,\mathbf{k}})-f(\epsilon_{j,\mathbf{k}})]\delta(\epsilon_{j,\mathbf{k}}-\epsilon_{i,\mathbf{k}}-\hbar\omega)
//kappa=(L22-L12*L21/L11)/eT^2,sigma=L11,
void writesigma(double*,double*,double*,int,int,int,int);
static double sum_factor;
double sigmod(const double& T, const double& mu, const double &E);

void Ele_Conductivity::Routine()
{
	if(INPUT.cond_method == 1)
	{
		this->method1();
	}
	else
	{
		this->method2();
	}
	return;
}

void Ele_Conductivity::method1()
{
	//get current work directory 
	char *wdbuf;
	wdbuf=getcwd(nullptr,0);
	//cout<<wdbuf<<endl;

	//Init
	const double pi=M_PI;
	int nf=INPUT.n_fwhm;//Calculate different fwhm at the same time.
	double st[nf], eps[nf];
	if(INPUT.smear==0)
		nf=1;
	else if(INPUT.smear==1)
	{
		for(int i=0;i<nf;i++)
			st[i]=INPUT.fwhm[i]/TWOSQRT2LN2;
	}
	else if(INPUT.smear==2)
	{
		for(int i=0;i<nf;i++)
			eps[i]=INPUT.fwhm[i]/2;
	}
	double factor = 2*pi/3*pow(P_HBAR,3)/(pow(P_ME,2)*pow(P_BOHR*1e-10,5));
	int allcount=-1;

	int ns=INPUT.nscf;//decide dynamic space of sigma_all, L22_all and L12_all
	int nfolder=INPUT.nscf;//number of folders containing different scf results.
	bool multi=true;
	double tot_sum=0;
	if(ns<=0)
	{ 
		ns=nfolder=1;
		multi=false;
	}
	if(!INPUT.error_con) ns=1;
	int nw=INPUT.wcut/INPUT.dw+1;
	double *sigma_all,*L22_all,*L12_all;
	int nsnfnw=nf*ns*nw;
	int nfnw=nf*nw;
	L22_all=new double[nsnfnw];
	L12_all=new double[nsnfnw];
	sigma_all=new double[nsnfnw];
	for(int i=0;i<nsnfnw;i++)
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
			std::cout<<"No "<<INPUT.multi_directory<<" exists."<<std::endl;
			exit(0);
		}
		chdir(INPUT.multi_directory.c_str());
	}
	
	//loop of snapshot(scf)
	for(int ifolder=0;ifolder<nfolder;ifolder++)
	{
	int is=ifolder%ns;
	if(multi)
	{
		int ifoldername=ifolder+1;
		string ifolderstr=int2str(ifoldername);
		string command="test -e "+ifolderstr;
		int exist=system(command.c_str());
		if(exist!=0)
		{
			std::cout<<"No "<<ifolderstr<<" exists."<<std::endl;
			exit(0);
		}
		chdir(ifolderstr.c_str());
	}

	Wavefunc WF;
	WfRead wfr(WF);
	wfr.Init();
	sum_factor=2*P_ME*INPUT.vol*pow(P_BOHR,3)/M_PI/P_QE/INPUT.nele/1e30/P_HBAR;
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
		std::cout<<"scf "<<ifolder<<" ; kpoint "<<ik+1<<std::endl;
		//WF.print(0);

		int nbb = (nband-1) * nband / 2;
		if(INPUT.cond_intra) nbb += nband;
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
		std::cout<<"nband: "<<nband<<std::endl;
			

		double factor2=factor/INPUT.dw*pow(WF.factor,2)/INPUT.vol;
		
		//loop of band
#ifdef _OPENMP
#pragma omp parallel for reduction(+:tot_sum, sigma_all[:nsnfnw], L12_all[:nsnfnw], L22_all[:nsnfnw]) schedule(dynamic)
#endif	
		for(int ib=0;ib<nband;++ib)
		{
			for(int jb=ib;jb<nband;++jb)
	    	{
				int ijb;
				double energyj=WF.eigE[jb];
				double energyi=WF.eigE[ib];
				double dE=energyj-energyi;
				double docc=WF.occ[ib]-WF.occ[jb];
				if(!INPUT.cond_intra)
				{
					if(docc<=1e-10 || dE < 1e-5) continue;
					ijb = nband*ib - ib*(ib+1)/2 + (jb - ib - 1);
				}
				else
				{
					ijb = ib*nband - (ib + 1) * ib / 2 + jb;
				}
				assert(dE>=0);
				int iw=int(dE/INPUT.dw);
				if(iw>=nw) continue;

				double corr2 = vmatrix[ijb];
				double L11fact;
				if(dE < 1e-5)
				{
					double avg_en = 0.5*(energyj+energyi);
					double occij = sigmod(INPUT.temperature/P_EV2K, INPUT.fermiE, avg_en);
					L11fact=factor2*corr2*(occij - 1)*occij/INPUT.temperature*P_EV2K*WF.wk;
				}
				else
				{
					if(INPUT.smearinvw)
					{
						L11fact=factor2/dE*docc*corr2;
					}
					else
					{
						L11fact=factor2*docc*corr2;
					}
				}
				if(ib == jb) L11fact *= 0.5;

				double Eij_F=(energyj+energyi)/2-INPUT.fermiE;
				double L12fact=Eij_F*L11fact;
				double L22fact=pow(Eij_F,2)*L11fact;
				if(INPUT.smear==0)
                {
                    int isiw=is*nw+iw;
                    sigma_all[isiw]=sigma_all[isiw]+L11fact;
                    L12_all[isiw]=L12_all[isiw]-L12fact;
                    L22_all[isiw]=L22_all[isiw]+L22fact;
                }
				else if(INPUT.smear==1||INPUT.smear==2)
				{
				  	double smearf;
					double smearpre;
				  	
				  	for(int ifi=0;ifi<nf;ifi++)
				  	{
						double pre;
						double st2;
						double eps2;
						if(INPUT.smear==1)
						{
							st2=-0.5/pow(st[ifi],2);
							pre=1/sqrt(2*pi)/st[ifi]*INPUT.dw;
						}
						else
						{
							eps2=pow(eps[ifi],2);
							pre=1/pi*eps[ifi]*INPUT.dw;
						}

						for(int iv=0;iv<nw;iv++)
						{
							if(INPUT.smear==1)
							{
								if(INPUT.smearinvw || dE < 1e-5)
								{
									smearf=exp(pow(INPUT.dw*(iv+0.5)-dE,2)*st2) + exp(pow(dE+INPUT.dw*(iv+0.5),2)*st2);
									smearpre = pre * smearf;
								}
								else
								{
									smearf=exp(pow(INPUT.dw*(iv+0.5)-dE,2)*st2) - exp(pow(dE+INPUT.dw*(iv+0.5),2)*st2);
									smearpre = pre/(iv+0.5)/INPUT.dw * smearf;
								}
							}
							else
							{
								if(INPUT.smearinvw  || dE < 1e-5)
								{
									smearf = 1/(eps2+pow(INPUT.dw*(iv+0.5)-dE,2)) + 1/(eps2+pow(dE+INPUT.dw*(iv+0.5),2));
									smearpre = pre * smearf;
								}
								else
								{
									smearf = 1/(eps2+pow(INPUT.dw*(iv+0.5)-dE,2)) - 1/(eps2+pow(dE+INPUT.dw*(iv+0.5),2));
									smearpre = pre/(iv+0.5)/INPUT.dw * smearf;
								}
								
							}
							if(iv<iw)
							{
								if(smearf<1e-6) continue;
							}
							else
							{
								if(smearf<1e-6) break;
							}
							int isifiv=is*nfnw+ifi*nw+iv;
							sigma_all[isifiv]+=smearpre*L11fact;
							L12_all[isifiv]-=smearpre*L12fact;
							L22_all[isifiv]+=smearpre*L22fact;
						}
				  	}
				}
				if(jb != ib)
				{
					if(dE < 1e-5)
					{
						double avg_en = 0.5*(energyj+energyi);
						double occij = sigmod(INPUT.temperature/P_EV2K, INPUT.fermiE, avg_en);
						tot_sum += (occij - 1)*occij/INPUT.temperature*P_EV2K*(corr2*WF.factor);
					}
					else
					{
						tot_sum += docc*(corr2*WF.factor/dE);
					}
				}
				//cout<<docc<<' '<<w<<' '<<corr2<<' '<<INPUT.dw<<' '<<WF.factor<<' '<<INPUT.vol<<' '<<factor<<endl;
	    	}
		}
		wfr.clean();
		delete[] vmatrix;
	}
	
	wfr.cleanclass();
	
	if(multi)
	{
		chdir("../");
	}
	
	}

	if(multi)
	{
		chdir(wdbuf);
	}
	double tot_sum_fac=double(4)/3/P_ME/P_QE*P_HBAR*P_HBAR/P_BOHR/P_BOHR*1e20/nfolder/INPUT.nele;
	tot_sum*=tot_sum_fac;

#ifdef __MPI
	double * f_sigma_all, * f_L12_all, * f_L22_all;
	double f_tot_sum;
	if(RANK==0)
	{
		f_sigma_all=new double[nsnfnw];
		f_L12_all=new double[nsnfnw];
		f_L22_all=new double[nsnfnw];
	}
	MPI_Reduce(sigma_all,f_sigma_all,nsnfnw,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(L12_all,f_L12_all,nsnfnw,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(L22_all,f_L22_all,nsnfnw,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	MPI_Reduce(&tot_sum,&f_tot_sum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	
	if(RANK==0)
	{
		cout<<"Total sum rule(from sum formula): "<<f_tot_sum<<endl;
		writesigma(f_sigma_all,f_L12_all,f_L22_all,nf,nw,ns,nfolder);
	}
#else
	cout<<"Total sum rule(from sum formula): "<<tot_sum<<endl;
	writesigma(sigma_all,L12_all,L22_all,nf,nw,ns,nfolder);
#endif

	delete []sigma_all;
	delete []L12_all;
	delete []L22_all;
#ifdef __MPI
	if(RANK==0)
	{
		delete[]f_sigma_all;
		delete[]f_L12_all;
		delete[]f_L22_all;
	}
#endif
    free(wdbuf);
	return;
}

void writesigma(double *sigma_all,double*L12_all,double* L22_all,int nf,int nw,int ns,int nfolder)
{
	const double pi=M_PI;
	double T=INPUT.temperature;
	string lag;
	ofstream ofs,ofs2,ofs3;
	double nfnw=nf*nw;
	for(int ifi=0;ifi<nf;ifi++)
	{	
	double* sigma=new double[nw];
	double* std_sigma=new double[nw];
	double* kappa=new double[nw];
	double* std_kappa=new double[nw];
	double* L12=new double [nw];
	double* L22=new double [nw];
	double* avg_kappa = new double[nw];
	if(INPUT.smear==0) lag=dou2str(INPUT.fwhm[ifi])+"I";
	if(INPUT.smear==1) lag=dou2str(INPUT.fwhm[ifi])+"G";
	if(INPUT.smear==2) lag=dou2str(INPUT.fwhm[ifi])+"L";
	string sigmaname="sigma"+lag+".txt";
	string kappaname="kappa"+lag+".txt";
	string Onsagername="Onsager"+lag+".txt";
	for(int iw=0;iw<nw;iw++)
	{
		sigma[iw]=0;
		kappa[iw]=0;
		L12[iw]=0;
		L22[iw]=0;
		for(int is=0;is<ns;is++)
		{
			int isifiw=is*nfnw+ifi*nw+iw;
			sigma[iw]+=sigma_all[isifiw];
			kappa[iw]+=(L22_all[isifiw]-(L12_all[isifiw]*L12_all[isifiw])/sigma_all[isifiw]);
			L12[iw]+=L12_all[isifiw];
			L22[iw]+=L22_all[isifiw];
		}
		sigma[iw]/=nfolder;
		kappa[iw]=kappa[iw]/nfolder/T;
		L12[iw]/=nfolder;
		L22[iw]/=nfolder;
		avg_kappa[iw] = (L22[iw] - L12[iw] * L12[iw] / sigma[iw])/T;
	}
	double sum=0;
	ofs.open(sigmaname.c_str());
	ofs2.open(kappaname.c_str());
	ofs3.open(Onsagername.c_str());
	if(ns>1)
	{
		ofs<<setw(8)<<"## w(eV) "<<setw(20)<<"sigma(Sm^-1)"<<setw(20)<<"error(Sm^-1)"<<endl;
		ofs2<<setw(8)<<"## w(eV) "<<setw(20)<<"kappa(W(mK)^-1)"<<setw(20)<<"error(W(mK)^-1)"<<endl;
	}
	else
	{
		ofs<<setw(8)<<"## w(eV) "<<setw(20)<<"sigma(Sm^-1)"<<endl;
		ofs2<<setw(8)<<"## w(eV) "<<setw(20)<<"kappa(W(mK)^-1)"<<endl;
	}
	ofs3 << setw(8) << "## w(eV) " << setw(20) << "sigma(Sm^-1)" << setw(20) << "kappa(W(mK)^-1)" << setw(20)
            << "L12/e(Am^-1)" << setw(20) << "L22/e^2(Wm^-1)" << endl;
	for(int iw=0;iw<nw;iw++)
	{
	  if(ns>1)
	  {
		std_kappa[iw]=std_sigma[iw]=0;
		for(int is=0;is<ns;is++)
		{
			int isifiw=is*nfnw+ifi*nw+iw;
			std_sigma[iw]+=pow(sigma[iw]-sigma_all[isifiw],2);
			std_kappa[iw]+=pow(kappa[iw]-(L22_all[isifiw]-(L12_all[isifiw]*L12_all[isifiw])/sigma_all[isifiw])/T,2);
		}
		std_sigma[iw]=std_sigma[iw]/ns/(ns-1);
		std_sigma[iw]=sqrt(std_sigma[iw])*INPUT.tpk;
		std_kappa[iw]=std_kappa[iw]/ns/(ns-1);
		std_kappa[iw]=sqrt(std_kappa[iw])*INPUT.tpk;
		ofs<<setw(8)<<(iw+0.5)*INPUT.dw<<setw(20)<<sigma[iw]<<setw(20)<<std_sigma[iw]<<endl;
		ofs2<<setw(8)<<(iw+0.5)*INPUT.dw<<setw(20)<<avg_kappa[iw]<<setw(20)<<std_kappa[iw]<<endl;
	  }
	  else
	  {
		ofs<<setw(8)<<(iw+0.5)*INPUT.dw<<setw(20)<<sigma[iw]<<endl;
		ofs2<<setw(8)<<(iw+0.5)*INPUT.dw<<setw(20)<<avg_kappa[iw]<<endl;
	  }
		ofs3<<setw(8)<<(iw+0.5)*INPUT.dw<<setw(20)<<sigma[iw]<<setw(20)<<avg_kappa[iw]<<setw(20)<<L12[iw]<<setw(20)<<L22[iw]<<endl;
		sum+=sigma[iw];
	}
	sum=sum*INPUT.dw*sum_factor;
	cout<<INPUT.fwhm[ifi]<<" Sum rule(from integral formula) : "<<sum<<endl;
	delete[]sigma;
	delete[]std_sigma;
	delete[]kappa;
	delete[]std_kappa;
	delete[]L22;
	delete[]L12;
	delete[]avg_kappa; 
	ofs.close();
	ofs2.close();
	ofs3.close();
	}
	return;
}

double sigmod(const double& T, const double& mu, const double &E)
{
	return 1.0/(1.0+exp((E-mu)/T));
}

	
