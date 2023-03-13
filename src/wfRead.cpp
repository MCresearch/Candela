#include "wfRead.h"
#include <stdio.h>
#include "matrixmultip.h"
//created by qianrui 2020-2-5

WfRead::WfRead(Wavefunc &wf)
{
	wfpt=&wf;
}

void WfRead::setval(Wavefunc &wf)
{
	wfpt=&wf;
}
void WfRead::ignore(int ik)
{
	if(INPUT.wf_in_type=="QE1")//QE1 for QE v5.4 ;QE2 for QE v6.4
	{
		//There is nothing to do when data of each kpoint are stored in separate files.
	}
	else if(INPUT.wf_in_type=="PWmat")
	{
		wfpwmat.readOCC(*wfpt,ik);
		wfpwmat.readGKK(*wfpt,ik);
		wfpwmat.readWF(*wfpt,ik);
	}
	else if(INPUT.wf_in_type=="ABACUS")
	{
		wfabacus.readOCC(*wfpt, ik);
	}
	else if(INPUT.wf_in_type=="QE2")
	{
		wfqe.readOCC2(*wfpt, ik);
	}
	else
	{
		cout<<"No such file type!"<<endl;
		exit(0);
	}
	return;
	
}

//read before the loop of ik
void WfRead:: Init()
{
	if(wfpt==nullptr)
	{
		cout<<"Class WfRead should be assigned."<<endl;
		exit(0);
	}
	if(INPUT.wf_in_type=="QE1")
	{
		wfqe.readOUT(*wfpt);
	}
	else if(INPUT.wf_in_type=="QE2")
	{
		wfqe.readOUT2(*wfpt);
	}
	else if(INPUT.wf_in_type=="PWmat")
	{
		wfpwmat.Init(*wfpt);	
	}
	else if(INPUT.wf_in_type=="ABACUS")
	{
		wfabacus.readOUT(*wfpt);
	}
	else
	{
		cout<<"No such file type: "<<INPUT.wf_in_type<<endl;
		exit(0);
	}
	return;
}

//read occ,gkk,wf under the loops of ik
void WfRead::readWF(int ik)
{
	
	if(wfpt==nullptr)
	{
		cout<<"Class WfRead should be assigned."<<endl;
		exit(0);
	}
	if(INPUT.wf_in_type=="QE1")
	{
		wfqe.readOCC(*wfpt,ik);
		wfqe.readGKK(*wfpt,ik);
		wfqe.readWF(*wfpt,ik);
	}
	else if(INPUT.wf_in_type=="QE2")
	{
		wfqe.readOCC2(*wfpt,ik);
		wfqe.readWF2(*wfpt,ik);
	}
	else if(INPUT.wf_in_type=="PWmat")
	{
		wfpwmat.readOCC(*wfpt,ik);
		wfpwmat.readGKK(*wfpt,ik);
		wfpwmat.readWF(*wfpt,ik);
	}
	else if(INPUT.wf_in_type=="ABACUS")
	{
		wfabacus.readOCC(*wfpt,ik);
		if(!INPUT.nonlocal)	wfabacus.readWF(*wfpt,ik);
	}
	else
	{
		cout<<"No such file type!"<<endl;
		exit(0);
	}
	return;
}

//read occ only
void WfRead::readOCC(int ik)
{
	
	if(wfpt==nullptr)
	{
		cout<<"Class WfRead should be assigned."<<endl;
		exit(0);
	}
	if(INPUT.wf_in_type=="QE1")
	{
		wfqe.readOCC(*wfpt,ik);
	}
	else if(INPUT.wf_in_type=="QE2")
	{
		wfqe.readOCC2(*wfpt,ik);
	}
	else if(INPUT.wf_in_type=="PWmat")
	{
		wfpwmat.readOCC(*wfpt,ik);
	}
	else if(INPUT.wf_in_type=="ABACUS")
	{
		wfabacus.readOCC(*wfpt,ik);
	}
	else
	{
		cout<<"No such file type!"<<endl;
		exit(0);
	}
	return;
}
void WfRead:: readvmatrix(const int ik, double* vmatrix)
{
	int nband = wfpt->nband;
	int nbb = (nband-1) * nband / 2;
	stringstream ss;
    ss<<INPUT.wfdirectory<<"vmatrix"<<ik+1<<".dat";
	binfstream binfs(ss.str(), "r");
	int head, tail;
	binfs>>head;
	ifnecheckv(nbb*8, head);
	rwread(binfs, vmatrix, nbb);
	binfs>>tail;
	ifnecheckv(head, tail);
}

void WfRead::calvmatrix(double* vmatrix)
{
	const int nband = wfpt->nband;
	const int nbb = (nband-1) * nband / 2;
	const int npw = wfpt->ngtot;
	ZEROS(vmatrix, nbb);
	complex<double> *pij = new complex<double>[nbb];
	complex<double> *pwave = new complex<double>[npw * nband];
	for (int id = 0; id < 3; ++id)
    {
		double *gid, kid;
		if(id == 0)
		{
			gid = wfpt->gkk_x;
			kid = wfpt->kpoint_x;
		}
		else if(id == 1)
		{
			gid = wfpt->gkk_y;
			kid = wfpt->kpoint_y;
		}
		else
		{
			gid = wfpt->gkk_z;
			kid = wfpt->kpoint_z;
		}

		// pxyz|right>
		for(int ib = 0 ; ib < nband ; ++ib)
		{
			for(int ig = 0 ; ig < npw ; ++ig)
			{
    	        pwave[ig + ib*npw] = wfpt->Wavegg[ib*npw + ig] * (gid[ig]+kid);
			}
		}
		dtrimultipAHB(nband,nband,npw, wfpt->Wavegg, npw, pwave, npw, pij, 1);
		if(INPUT.gamma)
		{
			for(int i = 0; i < nbb ; ++i)
			{
				vmatrix[i] += 4 * pow(pij[i].imag(), 2);
			}
		}
		else
		{
			for(int i = 0; i < nbb ; ++i)
			{
				vmatrix[i] += norm(pij[i]);
			}
		}
	}
	
	delete[] pij;
	delete[] pwave;
}

//act before the end of each loop
void WfRead:: clean()
{
	pclean(wfpt->gkk_x);
	pclean(wfpt->gkk_y);
	pclean(wfpt->gkk_z);
	pclean(wfpt->Wavegg);
	pclean(wfpt->occ);
	pclean(wfpt->eigE);
}
//act after the loop of ik
void WfRead::cleanclass()
{
	wfpwmat.clean();
	wfqe.clean();
}



