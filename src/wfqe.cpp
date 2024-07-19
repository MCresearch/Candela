#include "wfqe.h"
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
	
void WfQE::clean()
{
	pclean(wk);
}

void WfQE::readOUT2(Wavefunc & wf) 
{
	//open outfile
	string outname=INPUT.wfdirectory+"/data-file-schema.xml";
	ifskwt.open(outname.c_str());
	if(!ifskwt)
	{
		cout<<"Erro in opening OUT file!"<<endl;
		exit(0);
	}
	//cout<<outname<<" has been opened."<<endl;
	string useless;
	string txt;
	
	//get lattice parameter
	searchead(ifskwt,txt,"<atomic_structure",1);
	string alatstr=findstr(txt,"alat");
	alat=str2dou(alatstr);//in bohr
	//cout<<"lattice parameter: "<<alat<<" a.u."<<endl;

	//get crystal axes
	searchead(ifskwt,txt,"<cell>",1);
	for(int i=0;i<3;i++)
	{
		ifskwt>>useless;
		size_t st=useless.find('>',0);
		size_t le=useless.size()-st-1;
		string vxstr=useless.substr(st+1,le);
		axes[i].x=str2dou(vxstr);
		ifskwt>>axes[i].y;
		ifskwt>>useless;
		le=useless.find('<',0);
		string vzstr=useless.substr(0,le);
		axes[i].z=str2dou(vzstr);
		//cout<<axes[i].x<<' '<<axes[i].y<<' '<<axes[i].z<<endl;
	}
	INPUT.celldm1=axes[0].norm();
	INPUT.celldm2=axes[1].norm();
	INPUT.celldm3=axes[2].norm();
	double vol;
	vol=axes[0].x*(axes[1].y*axes[2].z-axes[2].y*axes[1].z)+ axes[0].y*(axes[2].x*axes[1].z-axes[1].x*axes[2].z)+ axes[0].z*(axes[1].x*axes[2].y-axes[2].x*axes[1].y);//vol in bohr^3
    vol=abs(vol);
    INPUT.vol=vol;
	
	//get gamma_only
	searchead(ifskwt,txt,"<basis>",1);
	ifskwt>>useless;
	if((useless=="<gamma_only>true</gamma_only>")||(useless=="<gamma_only>T</gamma_only>"))
		INPUT.gamma=true;
	else if(useless=="<gamma_only>false</gamma_only>"||useless=="<gamma_only>F</gamma_only>")
		INPUT.gamma=false;
	else
	{
		cout<<"Error in read GAMMA"<<endl;
		exit(0);
	}

	//get nband &nele
	searchead(ifskwt,txt,"<band_structure>",1);
	for(int i=0;i<4;i++)
		getline(ifskwt,txt);
	string nbandstr=findstr(txt);
	wf.nband=nband=str2int(nbandstr);
	getline(ifskwt,txt);	
	string nelestr=findstr(txt);
	INPUT.nele=str2dou(nelestr);
	//cout<<"number of electrons: "<<INPUT.nele<<endl;
	
	//get fermi energy
	for(int i=0;i<3;i++)
		getline(ifskwt,txt);
	string fermiEstr=findstr(txt);
	INPUT.fermiE=str2dou(nelestr);
	INPUT.fermiE*=P_HA;
	//get nkpoint
	searchead(ifskwt,txt,"</starting_k_points>",1);
	getline(ifskwt,txt);
	string nkstr=findstr(txt);
	INPUT.nkpoint=nkpoint=str2int(nkstr);
	
	//find the position of occ and eigen value
	searchead(ifskwt,txt,"<ks_energies>",1);
	wf.factor=1;
		
}

void WfQE::readOUT(Wavefunc & wf) 
{
	//open outfile
	string outname=INPUT.wfdirectory+"/data-file.xml";
	ifskwt.open(outname.c_str());
	if(!ifskwt)
	{
		cout<<"Erro in opening OUT file!"<<endl;
		exit(0);
	}
	//cout<<outname<<" has been opened."<<endl;
	string useless;
	string txt;

	//get lattice parameter
	searchead(ifskwt,txt,"<LATTICE_PARAMETER",1);
	ifskwt>>alat;//in bohr
	//cout<<"lattice parameter: "<<alat<<" a.u."<<endl;
	
	//get crystal axes
	searchead(ifskwt,txt,"<UNITS_FOR_DIRECT_LATTICE_VECTORS",1);
	getline(ifskwt,txt);
	ifskwt>>axes[0].x>>axes[0].y>>axes[0].z;
	for(int i=0;i<3;i++)
	{
		getline(ifskwt,txt);
	}
	ifskwt>>axes[1].x>>axes[1].y>>axes[1].z;
	for(int i=0;i<3;i++)
	{
		getline(ifskwt,txt);
	}
	ifskwt>>axes[2].x>>axes[2].y>>axes[2].z;
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
	searchead(ifskwt,txt,"<GAMMA_ONLY",1);
	string gamma;
	ifskwt>>gamma;
	if(gamma=="T") 
		INPUT.gamma=true;
	else if (gamma=="F")
		INPUT.gamma=false;
	else
	{
		cout<<"Error in read GAMMA"<<endl;
		exit(0);
	}

	//get weight of kpoint
	searchead(ifskwt,txt,"<NUMBER_OF_K-POINTS",1);
	ifskwt>>nkpoint;
	INPUT.nkpoint=nkpoint;
	for(int i=0;i<5;i++) getline(ifskwt,txt);
	wk=new double [nkpoint];
	for(int i=0;i<nkpoint;i++)
	{
		getline(ifskwt,txt);
		string wkstr=findstr(txt,"WEIGHT");
		wk[i]=str2dou(wkstr);
	}
	
	//get number of ele
	searchead(ifskwt,txt,"<NUMBER_OF_ELECTRONS",1);
	ifskwt>>INPUT.nele;
	//cout<<"number of electrons: "<<INPUT.nele<<endl;

	//get number of fermi energy
	searchead(ifskwt,txt,"<FERMI_ENERGY",1);
	ifskwt>>INPUT.fermiE;
	INPUT.fermiE*=P_HA;
	wf.factor=1;
	ifskwt.close();
}

void WfQE::readOCC(Wavefunc & wf, int & ik) 
{
	//open occfile
	string idwf=getidwf(ik);
	string occname=INPUT.wfdirectory+"/"+idwf+"/eigenval.xml";
	ifsocc.open(occname.c_str());
	if(!ifsocc)
	{
		cout<<"Erro in opening OCC file!"<<endl;
		exit(0);
	}
	//cout<<occname<<" has been opened."<<endl;
	string useless;
	string txt;

	//get occ & eigen Energy
	for(int i=0;i<9;i++)
		getline(ifsocc,txt);
	string nbndstr=findstr(txt,"size");
	nband=str2int(nbndstr);
	wf.nband=nband;
	wf.occ=new double [nband];
	wf.eigE=new double [nband];
	for(int i=0;i<nband;i++)
	{
		ifsocc>>wf.eigE[i];
		wf.eigE[i]*=P_HA;
	}
	for(int i=0;i<3;i++)
	{
		getline(ifsocc,txt);
	}
	wf.wk=wk[ik]/2;
	for(int i=0;i<nband;i++)
	{
		ifsocc>>wf.occ[i];
		wf.occ[i]=wf.occ[i]*wk[ik]/2;
	}
	ifsocc.close();
}

void WfQE::readOCC2(Wavefunc & wf, int & ik) 
{
	string useless;
	string txt;

	//get occ & eigen Energy
	getline(ifskwt,txt);
	string wkstr=findstr(txt,"weight");
	double kweight;
	kweight=str2dou(wkstr);
	getline(ifskwt,txt);
	useless=findstr(txt);
	ngtot=str2int(useless);
	wf.ngtot=ngtot;
	getline(ifskwt,txt);
	wf.occ=new double [nband];
	wf.eigE=new double [nband];
	for(int i=0;i<nband;i++)
	{
		ifskwt>>wf.eigE[i];
		wf.eigE[i]*=P_HA;
	}
	for(int i=0;i<3;i++)
	{
		getline(ifskwt,txt);
	}
	wf.wk=kweight/2;
	for(int i=0;i<nband;i++)
	{
		ifskwt>>wf.occ[i];
		wf.occ[i]=wf.occ[i]*kweight/2;
	}
	for(int i=0;i<4;i++)
	{
		getline(ifskwt,txt);
	}
}

void WfQE:: readGKK(Wavefunc & wf, int & ik)
{
	//open gkkfile and initialize
	string idwf=getidwf(ik);
	string gkkname=INPUT.wfdirectory+"/"+idwf+"/gkvectors.dat";
	rwsgkk.open(gkkname,"r");
	if(!rwsgkk)
	{
		cout<<"Error in opening GKK file!"<<endl;
		exit(0);
	}
	//cout<<gkkname<<" has been opened."<<endl;
	int strw,endrw,inttmp;
	string useless;
	locate(rwsgkk,useless,"NUMBER_OF_GK-VECTORS",1);
	
	//get ng (max num of PW)
	rwsgkk>>strw>>inttmp>>ngtot>>endrw;
	//cout<<"ngtot: "<<ngtot<<endl;
	wf.ngtot=ngtot;
	ifnecheckv(strw,endrw);
	locate(rwsgkk,useless,"K-POINT_COORDS",1);

	//get kpoint vector
	double kx_cry,ky_cry,kz_cry;
	rwsgkk>>strw>>inttmp>>kx_cry>>ky_cry>>kz_cry>>endrw;	
	wf.kpoint_x=kx_cry*2*M_PI/INPUT.celldm1;
	wf.kpoint_y=ky_cry*2*M_PI/INPUT.celldm2;
	wf.kpoint_z=kz_cry*2*M_PI/INPUT.celldm3;
	//cout<<"kpoint_vector: ("<<wf.kpoint_x<<','<<wf.kpoint_y<<','<<wf.kpoint_z<<")\n";
	ifnecheckv(strw,endrw);
	locate(rwsgkk,useless,"INDEX",1);
	int *idk,*ik_x,*ik_y,*ik_z;
	idk=new int [ngtot];
	ik_x=new int [ngtot];
	ik_y=new int [ngtot];
	ik_z=new int [ngtot];
	wf.gkk_x=new double [ngtot];
	wf.gkk_y=new double [ngtot];
	wf.gkk_z=new double [ngtot];

	//get index of each plane wave but temporarily unused
	rwsgkk>>strw>>inttmp;
	rwread(rwsgkk,idk,ngtot);
	delete[] idk;//Temporarily unused
	rwsgkk>>endrw;
	ifnecheckv(strw,endrw);
	locate(rwsgkk,useless,"GRID",1);

	//get plane wave
	rwsgkk>>strw>>inttmp;
	inttmp=((ngtot)*3+1)*4;
	ifnecheckv(strw,inttmp);
	for(int i=0;i<ngtot;i++)
	{
		rwsgkk>>ik_x[i]>>ik_y[i]>>ik_z[i];
		wf.gkk_x[i]=ik_x[i]*2*M_PI/INPUT.celldm1;
		wf.gkk_y[i]=ik_y[i]*2*M_PI/INPUT.celldm2;
		wf.gkk_z[i]=ik_z[i]*2*M_PI/INPUT.celldm3;
	}
	rwsgkk>>endrw;
	ifnecheckv(strw,endrw);
	delete []ik_x;
	delete []ik_y;
	delete []ik_z;
}
void WfQE::readWF(Wavefunc &wf, int &ik)
{
	//open wffile
	string idwf=getidwf(ik);
	string wfname=INPUT.wfdirectory+"/"+idwf+"/evc.dat";
	rwswf.open(wfname,"r");
	if(!rwswf)
	{
		cout<<"Erro in opening WF file!"<<endl;
		exit(0);
	}
	//cout<<wfname<<" has been opened."<<endl;
	int strw,endrw,inttmp;
	string useless;
	//get nband, nk, etc.
	locate(rwswf,useless,"INFO",1);
	
	string igwxstr=findstr(useless,"igwx");
	string nbndstr=findstr(useless,"nbnd");
	string ikstr=findstr(useless,"ik");
	string nkstr=findstr(useless,"nk");
	int ngtot_2=str2int(igwxstr);
	ifnecheckv(ngtot,ngtot_2);
	int nband_2=str2int(nbndstr);
	ifnecheckv(nband,nband_2);
	int ik_2=str2int(ikstr)-1;
	ifnecheckv(ik,ik_2);
	int nkpoint_2=str2int(nkstr);
	ifnecheckv(nkpoint,nkpoint_2);
	wf.Wavegg=new complex<double>[nband*ngtot];
	//double *sum;
	//sum=new double [nband];
	
	//read WF	
	for(int i=0;i<nband;i++)
	{
		//sum[i]=0;
		locate(rwswf,useless,"evc",1);
		rwswf>>strw>>inttmp;
		for(int j=0;j<ngtot;j++)
		{
			int index=i*ngtot+j;
			rwswf>>wf.Wavegg[index];
			//sum[i]+=pow(wf.Wavegg[index].real(),2)+pow(wf.Wavegg[index].imag(),2);
		}
		rwswf>>endrw;
		ifnecheckv(strw,endrw);
		locate(rwswf,useless,"evc",1);
	}
	for(int i=0;i<nband;i++)
	{
		//cout<<"iband: "<<i<<"\tsum : "<<sum[i]<<endl;
		//cout<<"iband "<<i<<" read"<<endl;
		/*for(int j=0;j<ngtot;j++)
			cout<<wf.Wavegg[i*ngtot+j]<<' ';
		cout<<endl;*/
	}
	
}

void WfQE::readWF2(Wavefunc &wf, int &ik)
{
	//open wffile
	int ikk=ik+1;
	string idwf=int2str(ikk);
	string wfname=INPUT.wfdirectory+"/wfc"+idwf+".dat";
	rwswf.open(wfname,"r");
	if(!rwswf)
	{
		cout<<"Erro in opening WF file!"<<endl;
		exit(0);
	}
	int strw,endrw,inttmp;
	int ik_2, ispin, gammaonly; 
	double scale;
	rwswf>>strw;
	//ik kx ky kz ispin gammaonly scale
	rwswf>>ik_2>>wf.kpoint_x>>wf.kpoint_y>>wf.kpoint_z>>ispin>>gammaonly>>scale;
	rwswf>>endrw;
	ifnecheckv(ik_2,ikk);
	ifnecheckv(strw,endrw);

	int npwx, ngtot_2, npol, nband_2;
	//npwx npw npol nbnd
	rwswf>>strw>>npwx>>ngtot_2>>npol>>nband_2>>endrw;
	ifnecheckv(strw,endrw);
	ifnecheckv(ngtot,ngtot_2);
	ifnecheckv(nband,nband_2);
	double b[9];

	rwswf>>strw;
	rwread(rwswf, b, 9);
	rwswf>>endrw;

	wf.gkk_x=new double [ngtot];
	wf.gkk_y=new double [ngtot];
	wf.gkk_z=new double [ngtot];
	wf.Wavegg=new complex<double> [nband*ngtot];
	//get plane wave
	rwswf>>strw;
	for(int i=0;i<ngtot;i++)
	{
		int gx,gy,gz;
		rwswf>>gx>>gy>>gz;
		wf.gkk_x[i] = gx*b[0] + gy*b[3] + gz*b[6];
		wf.gkk_y[i] = gx*b[1] + gy*b[4] + gz*b[7];
		wf.gkk_z[i] = gx*b[2] + gy*b[5] + gz*b[8];
		//cout<<wf.gkk_x[i]<<' '<<wf.gkk_y[i]<<' '<<wf.gkk_z[i]<<endl;
	}
	rwswf>>endrw;
	ifnecheckv(strw,endrw);
	for(int i=0;i<nband;i++)
	{
		rwswf>>strw;
		for(int j=0;j<ngtot;j++)
		{
			rwswf>>wf.Wavegg[i*ngtot+j];
		}
		rwswf>>endrw;
		ifnecheckv(strw,endrw);
	}
	
	rwswf.close();
	return;
}


string findstr(string longstr,const string &shortstr)
{
	size_t pos=longstr.find(shortstr,0);
	size_t quol,quor;
	if(pos==string::npos)
	{
		cout<<"No specific string is found."<<endl;
		exit(0);
	}
	else
	{
		quol=longstr.find('\"',pos+1);
		quor=longstr.find('\"',quol+1);
		if(quol==string::npos||quor==string::npos)
		{
			cout<<"No specific string is found."<<endl;
			exit(0);
		}
	}
	size_t length=quor-quol-1;
	return longstr.substr(quol+1,length);
} 

string findstr(string longstr)
{
	size_t quol,quor;
	quol=longstr.find('>',0);
	quor=longstr.find('<',quol+1);
	if(quol==string::npos||quor==string::npos)
	{
		cout<<"No specific string is found."<<endl;
		exit(0);
	}
	size_t length=quor-quol-1;
	return longstr.substr(quol+1,length);
}
void locate(binfstream &file,string& useless,const string obj,int n)
{
	int find=0;
	int inttmp;
	char str[1000];
	while(fgets(str,1000,file.fileptr))
	{
		useless.assign(str);
		size_t l;
		l=useless.find(obj,0);
		if(l!=string::npos)
		{
			find++;
		}
		if(find>=n)
		{ 
			fread(&inttmp,4,1,file.fileptr);
			break;
		}
       // cout<<useless;//for test
	}
	if(feof(file.fileptr))
	{
		cout<<"Can't find  "<<obj<<". File isn't enough!"<<endl;
        exit(0);
	}
}
 
string getidwf(int &ik)
{
	int ii=ik+1;//ik starts from 0;
	if(ii<10)
	{
		string st="K";
		stringstream ss;
		ss<<st<<"0000"<<ii;
		string result;
		ss>>result;
		return result;
	}
	else if(ii<100)
	{
		string st="K";
		stringstream ss;
		ss<<st<<"000"<<ii;
		string result;
		ss>>result;
		return result;
	}
	else if(ii<1000)
	{
		string st="K";
		stringstream ss;
		ss<<st<<"00"<<ii;
		string result;
		ss>>result;
		return result;
	}	
	else if(ii<10000)
	{
		string st="K";
		stringstream ss;
		ss<<st<<"0"<<ii;
		string result;
		ss>>result;
		return result;
	}
	else if(ii<100000)
	{
		string st="K";
		stringstream ss;
		ss<<st<<ii;
		string result;
		ss>>result;
		return result;
	}
	else
	{
		cout<<"The number of kpoints should be no more than 100000!!"<<endl;
		exit(0);
	}
}


