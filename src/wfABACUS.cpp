#include "wfABACUS.h"
#include "const.h"
#include "input.h"
#include <stdio.h>
#include <string>
#include <string.h>
#include <sstream>

//created by qianrui on 2020-2-5

string getidwf( int &);
void locate(binfstream& ,string&,const string,int);
string findstr(string ,const string &); 
string findstr(string); 
	
void WfABACUS::clean()
{
	pclean(wk);
	pclean(occ);
	pclean(energy);
}


void WfABACUS::readOUT(Wavefunc & wf) 
{
	//open outfile
	string outname=INPUT.wfdirectory+"/running_scf.log";
	ifskwt.open(outname.c_str());
	if(!ifskwt)
	{
		cout<<"Erro in opening OUT file!"<<endl;
		exit(0);
	}
	//cout<<outname<<" has been opened."<<endl;
	string useless;
	string txt;
	string checkstr;

	//get lattice parameter
	searchead(ifskwt,txt,"ntype",1);
	ifskwt>>useless>>useless>>txt>>useless>>alat;//in bohr
	checkstr="(Bohr)";
	ifnecheckv(txt,checkstr);
	//cout<<"lattice parameter: "<<alat<<" a.u."<<endl;
	
	//get crystal axes
	searchead(ifskwt,txt,"Lattice",1);
	ifskwt>>axes[0].x>>axes[0].y>>axes[0].z;
	ifskwt>>axes[1].x>>axes[1].y>>axes[1].z;
	ifskwt>>axes[2].x>>axes[2].y>>axes[2].z;
	for(int i=0;i<3;i++)
	{
		axes[i]=alat*axes[i];
	}
	/*cout<<"axes vector:"<<endl;
	for(int i=0;i<3;i++)
    {
		cout<<axes[i].x<<'\t'<<axes[i].y<<'\t'<<axes[i].z<<endl;
	}*/
	INPUT.celldm1=axes[0].norm();
	INPUT.celldm2=axes[1].norm();
	INPUT.celldm3=axes[2].norm();
	double vol;
	vol=axes[0].x*(axes[1].y*axes[2].z-axes[2].y*axes[1].z)+ axes[0].y*(axes[2].x*axes[1].z-axes[1].x*axes[2].z)+ axes[0].z*(axes[1].x*axes[2].y-axes[2].x*axes[1].y);//vol in bohr^3
    vol=abs(vol);
    INPUT.vol=vol;
	
	//get GAMMA_ONLY
		INPUT.gamma=false;

	//get number of ele
	searchead(ifskwt,txt,"electron",1);
	ifskwt>>useless>>txt>>useless>>useless>>useless>>useless>>useless>>INPUT.nele;
	checkstr="electron";
	ifnecheckv(txt,checkstr);
	//cout<<"number of electrons: "<<INPUT.nele<<endl;

	//get number of kpoint
	searchead(ifskwt,txt,"nkstot_ibz",1);
	nkpoint = str2int(txt.substr(2,10));
	ifskwt >> txt;
	ifnecheckv(txt,string("IBZ"));
	
	searchead(ifskwt,txt,"KPOINTS",1);
	INPUT.nkpoint=nkpoint;
	wk=new double [nkpoint];

	for(int i=0;i<nkpoint;i++)
	{
		ifskwt>>useless>>useless>>useless>>useless>>wk[i];
		//cout<<"wk"<<wk[i]<<endl;
	}

	//get nband
	searchead(ifskwt,txt,"occupied",1);
	ifskwt>>txt>>useless>>nband;
	checkstr="NBANDS";
	ifnecheckv(checkstr,txt);
	wf.nband=nband;
	

	//get number of fermi energy
	searchead(ifskwt,txt,"E_exx",1);
	ifskwt>>txt>>useless>>INPUT.fermiE; //in eV
	checkstr="E_Fermi";
	ifnecheckv(checkstr,txt);
	wf.factor=1;
	ifskwt.close();
	
	//get occ
	energy=new double[nkpoint*nband];
	occ=new double [nkpoint*nband];
	string occname=INPUT.wfdirectory+"/istate.info";
	ifsocc.open(occname.c_str());
	if(!ifsocc)
	{
		cout<<"Erro in opening OCC file!"<<endl;
		exit(0);
	}
	for(int ik = 0; ik < nkpoint; ++ik)
	{
		ifsocc>>useless;
		getline(ifsocc,txt);
		for(int i=0;i<nband;i++)
		{
			ifsocc>>useless>>energy[ik*nband+i]>>occ[ik*nband+i];
		}
	}
	ifsocc.close();

}

void WfABACUS::readOCC(Wavefunc & wf, int & ik) 
{
	wf.occ=new double [nband];
	wf.eigE=new double [nband];
	wf.wk=wk[ik];
	for(int i=0;i<nband;i++)
	{
		wf.eigE[i]=energy[ik*nband+i];
		wf.occ[i]=occ[ik*nband+i]/2;
	}
}


void WfABACUS::readWF(Wavefunc &wf, int &ik)
{
	//open wavefunc file and initialize
	string wfname=INPUT.wfdirectory+"/WAVEFUNC"+int2str(ik+1)+".dat";
	rwswf.open(wfname,"r");
	if(!rwswf)
	{
		cout<<"Error in opening WAVEFUNC*.dat file!"<<endl;
		exit(0);
	}

	int strw,endrw;
	int ik_2,nkpoint_2,nband_2;
	double lat0,invlat0,wk_2,ecut;
	double kx_cry,ky_cry,kz_cry;
	
	//get data
	rwswf>>strw>>ik_2>>nkpoint_2>>kx_cry>>ky_cry>>kz_cry>>wk_2>>ngtot>>nband_2>>ecut>>lat0>>invlat0>>endrw;
	ik_2-=1;
	ifnecheckv(ik,ik_2);
	ifnecheckv(nkpoint_2,nkpoint);
	ifnecheckv(nband,nband_2);
	ifnecheckv(strw,endrw);
	//cout<<"ngtot: "<<ngtot<<endl;
	wf.ngtot=ngtot;

	//get kpoint vector	
	wf.kpoint_x=kx_cry*invlat0;
	wf.kpoint_y=ky_cry*invlat0;
	wf.kpoint_z=kz_cry*invlat0;
	//cout<<"kpoint_vector: ("<<wf.kpoint_x<<','<<wf.kpoint_y<<','<<wf.kpoint_z<<")\n";

	//get inverse lattice matrix
	double e11,e12,e13,e21,e22,e23,e31,e32,e33;
	rwswf>>strw>>e11>>e12>>e13>>e21>>e22>>e23>>e31>>e32>>e33>>endrw;
	ifnecheckv(strw,endrw);
	
	//get gkk
	rwswf>>strw;
	double *gcar_x,*gcar_y,*gcar_z;
	gcar_x=new double [ngtot];
	gcar_y=new double [ngtot];
	gcar_z=new double [ngtot];
	wf.gkk_x=new double [ngtot];
	wf.gkk_y=new double [ngtot];
	wf.gkk_z=new double [ngtot];

	int inttmp=((ngtot)*3)*8;
	ifnecheckv(strw,inttmp);

	wf.ig0 = -1;
	for(int i=0;i<ngtot;i++)
	{
		rwswf>>gcar_x[i]>>gcar_y[i]>>gcar_z[i];
		if(pow(gcar_x[i],2)+pow(gcar_y[i],2)+pow(gcar_z[i],2) < 1e-8)
		{
			wf.ig0 = i;
			// cout<<"ig0: "<<wf.ig0<<endl;
		}
		// cout<<gcar_x[i]<<' '<<gcar_y[i]<<' '<<gcar_z[i]<<endl;
		wf.gkk_x[i]=gcar_x[i]*invlat0;
		wf.gkk_y[i]=gcar_y[i]*invlat0;
		wf.gkk_z[i]=gcar_z[i]*invlat0;
	}
	rwswf>>endrw;
	ifnecheckv(strw,endrw);

	delete []gcar_x;
	delete []gcar_y;
	delete []gcar_z;

	//read wavefunc
	wf.Wavegg=new complex<double>[nband*ngtot];
	//double *sum;
	//sum=new double [nband];
	
	//read WF	
	for(int i=0;i<nband;i++)
	{
		//sum[i]=0;
		rwswf>>strw;
		for(int j=0;j<ngtot;j++)
		{
			int index=i*ngtot+j;
			rwswf>>wf.Wavegg[index];
			//sum[i]+=pow(wf.Wavegg[index].real(),2)+pow(wf.Wavegg[index].imag(),2);
		}
		rwswf>>endrw;
		ifnecheckv(strw,endrw);
	}
	for(int i=0;i<nband;i++)
	{
		//cout<<"iband: "<<i<<"\tsum : "<<sum[i]<<endl;
		//cout<<"iband "<<i<<" read"<<endl;
		/*for(int j=0;j<ngtot;j++)
			cout<<wf.Wavegg[i*ngtot+j]<<' ';
		cout<<endl;*/
	}
	return;
}



 


