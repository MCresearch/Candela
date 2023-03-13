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
double Velocity_Matrix_Local(const int,const int,Wavefunc &);
//double Velocity_Matrix_NonLocal(const int,const int,Wavefunc &,const int);
//double CPU_Velocity_Matrix_NonLocal(const int,const int,Wavefunc &);

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
	double w;
	int iw;
	complex<double> corrx,corry,corrz;
	double corr2;
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
	int is=ifolder%ns;
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
		cout<<"scf "<<ifolder<<" ; kpoint "<<ik+1<<endl;
		cout<<"nband: "<<nband<<endl;
		//WF.print(0);

		const int nbb = (nband-1) * nband / 2;
		double *vmatrix = new double [nbb];
		if(INPUT.nonlocal)
		{
			wfr.readvmatrix(ik, vmatrix);
		}
		else
		{
			for(int ib=0;ib<nband;ib++)
			{
				WF.checknorm(ik,ib);//check if wavefunction is normalized to 1.
			}
			wfr.calvmatrix(vmatrix);
		}
			

		double factor2=factor/INPUT.dw*pow(WF.factor,2)/INPUT.vol;
		
		//loop of band
		int ijb = 0;
		for(int ib=0;ib<nband;++ib)
		{
			for(int jb=ib+1;jb<nband;++jb, ++ijb)
	    	{
				double energyj=WF.eigE[jb];
				double energyi=WF.eigE[ib];
				w=energyj-energyi;
				double docc=WF.occ[ib]-WF.occ[jb];
				if(docc<=0.0000000001) continue;
				assert(w>=0);
				iw=int(w/INPUT.dw);
				if(iw>=nw||w==0) continue;//add w==0, sometimes w can be 0 due to truncation error.

				corr2 = vmatrix[ijb];
				//corr2=Velocity_Matrix_Local(ib,jb,WF); 
				
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
								if(INPUT.smearinvw)
								{
									smearf=exp(pow(INPUT.dw*(iv+0.5)-w,2)*st2)+exp(pow(w+INPUT.dw*(iv+0.5),2)*st2);
									smearpre = pre * smearf;
								}
								else
								{
									smearf=exp(pow(INPUT.dw*(iv+0.5)-w,2)*st2) - exp(pow(w+INPUT.dw*(iv+0.5),2)*st2);
									smearpre = pre/(iv+0.5)/INPUT.dw * smearf;
								}
								
							}
							else
							{
								smearf = 1/(eps2+pow(INPUT.dw*(iv+0.5)-w,2))+1/(eps2+pow(w+INPUT.dw*(iv+0.5),2));
								smearpre=pre*smearf;
								
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
	
				tot_sum=tot_sum+(WF.occ[ib]-WF.occ[jb])*(corr2*WF.factor/(WF.eigE[jb]-WF.eigE[ib]));
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
	}
	double sum=0;
	ofs.open(sigmaname.c_str());
	ofs2.open(kappaname.c_str());
	ofs3.open(Onsagername.c_str());
	if(ns>1)
	{
		ofs<<"#w\tsigma\terror"<<endl;
		ofs2<<"#w\tkappa\terror"<<endl;
	}
	else
	{
		ofs<<"#w\tsigma"<<endl;
		ofs2<<"#w\tkappa"<<endl;
	}
		ofs3<<"#w\tL11\tL12/e\tL22/e^2"<<endl;
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
		ofs<<(iw+0.5)*INPUT.dw<<'\t'<<sigma[iw]<<'\t'<<std_sigma[iw]<<endl;
		ofs2<<(iw+0.5)*INPUT.dw<<'\t'<<kappa[iw]<<'\t'<<std_kappa[iw]<<endl;
		}
	  else
	  {
		ofs<<(iw+0.5)*INPUT.dw<<'\t'<<sigma[iw]<<endl;
		ofs2<<(iw+0.5)*INPUT.dw<<'\t'<<kappa[iw]<<endl;
	  }
		ofs3<<(iw+0.5)*INPUT.dw<<'\t'<<sigma[iw]<<'\t'<<L12[iw]<<'\t'<<L22[iw]<<endl;
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
	ofs.close();
	ofs2.close();
	ofs3.close();
	}
	return;
}

double Velocity_Matrix_Local(const int ib,const int jb,Wavefunc & WF)
{
		double corr2;
		int n=WF.ngtot;
		Vector3<double> gplusk,corrim(0,0,0);
		complex<double> corrx,corry,corrz,C_ji_c,C_ii,C_jiii;
		double C_jiii_im,corrxim,corryim,corrzim;
		corrx=corry=corrz=0;
		for(int ig=0;ig<n;ig++)
        {
			C_ji_c=conjugate(WF.Wavegg[jb*n+ig]);
			C_ii=WF.Wavegg[ib*n+ig];
			gplusk.x=WF.gkk_x[ig]+WF.kpoint_x;
			gplusk.y=WF.gkk_y[ig]+WF.kpoint_y;
			gplusk.z=WF.gkk_z[ig]+WF.kpoint_z;
			if(INPUT.gamma)
			{
				C_jiii_im=C_ji_c.real()*C_ii.imag()+C_ji_c.imag()*C_ii.real();
				corrim+=C_jiii_im*gplusk;
			}
			else
			{
				C_jiii=C_ji_c*C_ii;
				corrx+=C_jiii*gplusk.x;
				corry+=C_jiii*gplusk.y;
				corrz+=C_jiii*gplusk.z;
			}
        }
			
        if(INPUT.gamma)
        {
            //Note that gkk_x[0]=gkk_y[0]=gkk_z[0]=0. Thus, we directly multiply corr by 2.
			corr2=pow(corrim.x,2)+pow(corrim.y,2)+pow(corrim.z,2);
			corr2*=4;
        }
		else
		{
			corr2=pow(corrx.real(),2)+pow(corrx.imag(),2)+pow(corry.real(),2)+pow(corry.imag(),2)+pow(corrz.real(),2)+pow(corrz.imag(),2);
        }
		return corr2;
}

	
