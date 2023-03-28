#include "wavefunc.h"
#include "input.h"
#include <iostream>
Wavefunc::Wavefunc(Wavefunc & wf)
{
	this->Wavegg=wf.Wavegg;
	this->gkk_x=wf.gkk_x;
	this->gkk_y=wf.gkk_y;
	this->gkk_z=wf.gkk_z;
	this->occ=wf.occ;
	this->eigE=wf.eigE;
	this->nband=wf.nband;
	this->ngtot=wf.ngtot;
}
Wavefunc:: ~Wavefunc()
{
	pclean(gkk_x);
    pclean(gkk_y);
    pclean(gkk_z);
    pclean(Wavegg);
    pclean(occ);
    pclean(eigE);
}
void Wavefunc:: checknorm(int ik,int iband)
{
	double corr=0;
	if(!INPUT.gamma)
	{
		for(int i=0;i<ngtot;i++)
		{
			if(i==ig0) continue;
			corr+=pow(Wavegg[i+iband*ngtot].real(),2)+pow(Wavegg[i+iband*ngtot].imag(),2);
		}
	}
	else
	{
		for(int i=0;i<ngtot;i++)
		{
			if(i==ig0)
			{ 
				corr+=pow(Wavegg[i+iband*ngtot].real(),2)+pow(Wavegg[i+iband*ngtot].imag(),2);
			}
			else	corr+=2*(pow(Wavegg[i+iband*ngtot].real(),2)+pow(Wavegg[i+iband*ngtot].imag(),2));
		}
	}
	corr *= factor;
	if(abs(1-corr)>1e-6)
	{
		cout<<"1-sum: "<<1-corr<<endl;
		cout<<"kpoint "<<ik<<",Band "<<iband<<" is not normalized!"<<endl;
		exit(0);
	}

}
void Wavefunc:: print(int iband)
{
	cout<<"Printing WF..."<<endl;
	cout<<"nband: "<<nband<<", ngtot: "<<ngtot<<endl;
	cout<<"eigenvalue, occ"<<endl;
	for(int i=0;i<nband;i++)
	{
		cout<<eigE[i]<<'\t'<<occ[i]<<endl;
	}
	double corr=0;
	if(!INPUT.gamma)
	{
		for(int i=0;i<ngtot;i++)
			corr+=pow(Wavegg[i+iband*ngtot].real(),2)+pow(Wavegg[i+iband*ngtot].imag(),2);
		cout<<"sum: "<<corr*factor<<endl;
	}
	else
	{
		for(int i=0;i<ngtot;i++)
		{
			if(i==ig0)
			{ 
				corr+=pow(Wavegg[i+iband*ngtot].real(),2)+pow(Wavegg[i+iband*ngtot].imag(),2);
			}
			else	corr+=2*(pow(Wavegg[i+iband*ngtot].real(),2)+pow(Wavegg[i+iband*ngtot].imag(),2));
		}
		cout<<"sum for gamma_only: "<<corr*factor<<endl;
	}
	cout<<"BZ k: "<<kpoint_x<<' '<<kpoint_y<<' '<<kpoint_z<<endl;
	cout<<"Coe,|C|^2,(kx,ky,kz),k^2/2"<<endl;
	double sqfact = sqrt(factor);
	for(int i=0;i<ngtot;i++)
	{
		double kx=gkk_x[i]+kpoint_x;
		double ky=gkk_y[i]+kpoint_y;
		double kz=gkk_z[i]+kpoint_z;
		double k2=pow(kx,2)+pow(ky,2)+pow(kz,2);
		k2/=2;
		cout<<Wavegg[i+iband*ngtot]*sqfact<<' '<<(pow(Wavegg[i+iband*ngtot].real(),2)+pow(Wavegg[i+iband*ngtot].imag(),2))*factor<<" ( "<<kx<<" , "<<ky<<" , "<<kz<<" ) "<<k2<<endl;
	}
	cout<<"End"<<endl;
	
}
