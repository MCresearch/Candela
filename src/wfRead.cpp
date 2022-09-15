#include "wfRead.h"
#include <stdio.h>
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
	else if(INPUT.wf_in_type=="PWmat"||INPUT.wf_in_type=="QE2"||INPUT.wf_in_type=="ABACUS")
	{
		this->readWF(ik);//It is needed when data of all kpoints are stored in a single file. 	
	}
	else
	{
		cout<<"No such file type!"<<endl;
		exit(0);
	}
	return;
	
}

void WfRead::ignoreOCC(int ik)
{
	if(INPUT.wf_in_type=="QE1")//QE1 for QE v5.4 ;QE2 for QE v6.4
	{
		//There is nothing to do when occupied numbers and eigen energies of each kpoint are stored in separate files.
	}
	else if(INPUT.wf_in_type=="PWmat"||INPUT.wf_in_type=="QE2"||INPUT.wf_in_type=="ABACUS")
	{
		this->readOCC(ik);	//It is needed when occupied numbers and eigen energies of all kpoints are sotred in a single file.
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
		wfabacus.readWF(*wfpt,ik);
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



