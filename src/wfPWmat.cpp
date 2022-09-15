#include "input.h"
#include "wfPWmat.h"
#include "gfun.h"
#include "const.h"
//qianrui
WfPWmat::WfPWmat()
{
	ngtotnod_9=nullptr;
    energy=nullptr;
	occ=nullptr;
	wk=nullptr;
    init=false;
}
WfPWmat::~WfPWmat()
{
	pclean(ngtotnod_9);
	pclean(energy);
	pclean(occ);
	pclean(wk);
}
void WfPWmat::clean()
{
	pclean(ngtotnod_9);
	pclean(energy);
	pclean(occ);
	pclean(wk);
}
void WfPWmat::Init(Wavefunc & wf)
{
//OUT.GKK
	int strw,endrw;
	int is_SO,islda;
	double Ecut;
	string gkkfile=INPUT.wfdirectory+"/OUT.GKK";
	rwsgkk.open(gkkfile.c_str(),"r");
	if(!rwsgkk)
    {
        cout<<"Erro in opening GKK file!"<<endl;
        exit(0);
    }
	//cout<<gkkfile<<" has been opened."<<endl;
	rwsgkk>>strw>>n1>>n2>>n3>>mg_nx>>nnodes>>nkpoint>>is_SO>>islda>>endrw;
	INPUT.nkpoint=nkpoint;
	//cout<<"There are "<<nkpoint<<" kpoints"<<endl;
	ifnecheckv(strw,endrw)
	rwsgkk>>strw>>Ecut>>endrw;
	ifnecheckv(strw,endrw)
	rwsgkk>>strw;
	for(int i=0;i<3;i++)
	{
		rwsgkk>>axes[i].x>>axes[i].y>>axes[i].z;//need test
		axes[i]/=P_BOHR;//axes in bohr
	}
	rwsgkk>>endrw;
	ifnecheckv(strw,endrw)
	int nprod=nkpoint*nnodes;
	int *ngtotnod_9_t;
	ngtotnod_9=new int[nprod];
	ngtotnod_9_t=new int[nprod];
	int nnodes_t;
	rwsgkk>>strw>>nnodes_t;
	ifnecheckv(nnodes,nnodes_t)
	rwread(rwsgkk,ngtotnod_9,nprod);
	rwsgkk>>endrw;
	ifnecheckv(strw,endrw)

//OUT.WG
	string wgfile=INPUT.wfdirectory+"/OUT.WG";
	rwswf.open(wgfile.c_str(),"r");
	if(!rwswf)
    {
        cout<<"Erro in opening WG file!"<<endl;
        exit(0);
    }
    //cout<<wgfile<<" has been opened."<<endl;
	int n1_t,n2_t,n3_t,mg_nx_t,nkpoint_t;
	int is_SO_2,islda_2;
	double Ecut_t;
	rwswf>>strw>>n1_t>>n2_t>>n3_t>>nband>>mg_nx_t>>nnodes_t>>nkpoint_t>>is_SO_2>>islda_2>>endrw;
	wf.nband=nband;
	ifnecheckv(strw,endrw)
	ifnecheckv(mg_nx_t,mg_nx)
	ifnecheckv(is_SO_2,is_SO)
	ifnecheckv(islda_2,islda)
	ifnecheckv(n1,n1_t)
	ifnecheckv(n2,n2_t)
	ifnecheckv(n3,n3_t)
	ifnecheckv(nkpoint_t,nkpoint)
	rwswf>>strw>>Ecut_t>>endrw;
	ifnecheckv(strw,endrw)
	double axes_t[9];
	rwswf>>strw>>axes_t>>endrw;
	ifnecheckv(strw,endrw)
	rwswf>>strw>>nnodes_t;
	rwread(rwswf,ngtotnod_9_t,nprod);
	rwswf>>endrw;
	ifnecheckv(nnodes_t,nnodes)
	ifnecheckv(strw,endrw)
	ifnecheckp(ngtotnod_9,ngtotnod_9_t,nprod)
	delete[]ngtotnod_9_t;
//OUT.FERMI
	string useless;
	ifstream ifsfermi("OUT.FERMI");
	ifsfermi >>useless>>useless>>useless>>INPUT.fermiE;
	ifsfermi.close();
//REPORT
	wk=new double [nkpoint];
	string txt;
    ifstream ifsrep("REPORT");
	searchead(ifsrep,txt,"IN.SOLVENT",1);
	ifsrep>>useless>>useless>>INPUT.nele;
	searchead(ifsrep,txt,"total",1);
	for(int ik=0;ik<nkpoint;ik++)
	{
		ifsrep>>useless>>useless>>useless>>wk[ik];
	}
	ifsrep.close();
//OUT.OCC
	string occfile=INPUT.wfdirectory+"/OUT.OCC";
	ifsocc.open(occfile.c_str());
	if(!ifsocc)
    {
        cout<<"Erro in opening OCC file!"<<endl;
        exit(0);
    }
    //cout<<occfile<<" has been opened."<<endl;
	kpoint_x=new double[nkpoint];
	kpoint_y=new double[nkpoint];
	kpoint_z=new double[nkpoint];
	wf.kpoint_x=wf.kpoint_y=wf.kpoint_z=0;//For PWmat gkk has store the information of kpoint.
	energy=new double[nkpoint*nband];
	occ=new double [nkpoint*nband];
	int iiband;
	for(int i=0;i<nkpoint;i++)
	{
		
		ifsocc>>useless>>useless>>kpoint_x[i]>>kpoint_y[i]>>kpoint_z[i];
		if(ifsocc.eof())
		{
			cout<<"Error in occ file. Not enough kpoints!"<<endl;
			exit(0);
		}
		getline(ifsocc,useless);
		getline(ifsocc,useless);
		for(int j=0;j<nband;j++)
		{
			ifsocc>>iiband>>energy[i*nband+j]>>occ[i*nband+j];
			occ[i*nband+j]/=2;//!!!! The max occupation number in PWmat is 2. 
			if(j==nband-1) 
			{
				ifnecheckv(iiband,nband);
			}
		}
		getline(ifsocc,useless);
	}
	ifsocc.close();
	delete []kpoint_x;
	delete []kpoint_y;
	delete []kpoint_z;
//calculate Volume
	vol=axes[0].x*(axes[1].y*axes[2].z-axes[2].y*axes[1].z)+ axes[0].y*(axes[2].x*axes[1].z-axes[1].x*axes[2].z)+ axes[0].z*(axes[1].x*axes[2].y-axes[2].x*axes[1].y);//vol in bohr^3
	vol=abs(vol);
	INPUT.vol=vol;
	wf.factor=vol;
	init=true;
	return;
}
void WfPWmat:: readOCC(Wavefunc & wf,int &ik)
{
	if(!init)
	{
		cout<<"Plz Init first!"<<endl;
		exit(0);
	}	
	wf.eigE=new double [nband];
	wf.occ=new double [nband];
	wf.wk=wk[ik];
	for(int i=0;i<nband;i++)
	{
		wf.eigE[i]=energy[ik*nband+i];
		wf.occ[i]=occ[ik*nband+i];
	}
}

void WfPWmat:: readGKK(Wavefunc & wf,int &ik)
{
	int strw,endrw;
	if(!init)
	{
		cout<<"Plz Init first!"<<endl;
		exit(0);
	}
	double *gkk_x,*gkk_y,*gkk_z;
	double *gkk_n_tmp,*gkk_n_xtmp,*gkk_n_ytmp,*gkk_n_ztmp;
    gkk_x=new double [mg_nx*nnodes];
    gkk_y=new double [mg_nx*nnodes];
    gkk_z=new double [mg_nx*nnodes];
    gkk_n_tmp=new double[mg_nx];
    gkk_n_xtmp=new double[mg_nx];
    gkk_n_ytmp=new double[mg_nx];
    gkk_n_ztmp=new double[mg_nx];
	ngtot=0;
	for(int j=0;j<nnodes;j++)
	{
		rwsgkk>>strw;
		rwread(rwsgkk,gkk_n_tmp,mg_nx);
		rwsgkk>>endrw;
		ifnecheckv(strw,endrw)
		rwsgkk>>strw;
		rwread(rwsgkk,gkk_n_xtmp,mg_nx);
		rwsgkk>>endrw;
		ifnecheckv(strw,endrw)
		rwsgkk>>strw;
		rwread(rwsgkk,gkk_n_ytmp,mg_nx);
		rwsgkk>>endrw;
		ifnecheckv(strw,endrw)
		rwsgkk>>strw;
		rwread(rwsgkk,gkk_n_ztmp,mg_nx);
		rwsgkk>>endrw;
		ifnecheckv(strw,endrw)
		//cout<<"ngtotnod  " <<ngtotnod_9[ikpt*nnodes+j]<<endl;//test
		for(int k=0;k<ngtotnod_9[ik*nnodes+j];k++)
		{
			gkk_x[ngtot+k]=gkk_n_xtmp[k];
			gkk_y[ngtot+k]=gkk_n_ytmp[k];
			gkk_z[ngtot+k]=gkk_n_ztmp[k];
			//cout<<k<<' '<<gkk_n_tmp[k]<<' '<<gkk_n_ztmp[k]<<' '<<gkk_n_ytmp[k]<<' '<<gkk_n_ztmp[k]<<endl;//test
		}
		ngtot=ngtot+ngtotnod_9[ik*nnodes+j];
	}
	wf.ngtot=ngtot;
	wf.gkk_x=new double [ngtot];
	wf.gkk_y=new double [ngtot];
	wf.gkk_z=new double [ngtot];
	for(int igg=0;igg<ngtot;igg++)
	{
		wf.gkk_x[igg]=gkk_x[igg];
		wf.gkk_y[igg]=gkk_y[igg];
		wf.gkk_z[igg]=gkk_z[igg];
	}
	delete[]gkk_x;
	delete[]gkk_y;
	delete[]gkk_z;
	delete[]gkk_n_tmp;
	delete[]gkk_n_xtmp;
	delete[]gkk_n_ytmp;
	delete[]gkk_n_ztmp;
	return;
}

 void WfPWmat::readWF(Wavefunc &wf,int &ik)
{
	int strw,endrw;
	if(!init)
	{
		cout<<"Plz Init first!"<<endl;
		exit(0);
	}
	//OUT.WG
	complex<double> *tmp_ug;
	tmp_ug=new complex<double>[mg_nx];
	wf.Wavegg=new complex<double>[nband*ngtot];
	for(int i=0;i<nband;i++)
	{
		rwswf>>strw;
		rwread(rwswf,tmp_ug,mg_nx);//tmp_ug in bohr^-3/2
		rwswf>>endrw;
		for(int igg=0;igg<ngtot;igg++)
		{
			wf.Wavegg[i*ngtot+igg]=tmp_ug[igg];
		}
		ifnecheckv(strw,endrw)
	}
	delete[] tmp_ug;
	return;
}
