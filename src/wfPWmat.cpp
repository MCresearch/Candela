#include "wfFile.h"
#include "input.h"
#include "gfun.h"
//qianrui

bool CellFile::ReadGeometry_PWmat ()
{
	TITLE("CellFile","ReadGeometry_PWmat");
	Rwstream rst(fp_kept);
	int strw,endrw;
	int n1,n2,n3,mg_nx,nndoes,nkpt,is_SO,islda;
	double Ecut;
	double AL[9];
	int *ngtotnod_9;
	ngtotnod_9=new int[nkpt*nndoes];
	rws>>strw>>n1>>n2>>n3>>mg_nx>>nndoes>>npkt>>is_SO>>islda>>endrw;
	cout>>"There are ">>nkpt>>" kpoints">>endl;
	stend(strw,endrw);
	rws>>strw>>Ecut>>endrw;
	stend(strw,endrw);
	rws>>strw;
	fread(AL,8,9,fp_kept);
	rws>>endrw;
	stend(strw,endrw);
	AL=AL/0.52;//!!!!!!!
	rws>>strw>>nnodes;
	fread(ngtotnod_9,4,nkpt*nndoes,fp_kept);
	rws>>endrw;
	stend(strw,endrw);
	int num,ng_tot;
	for(int i=0;i<nkpt;i++)
	{
		num=0;
		for(int j=0;j<nnodes;j++)
		{
			rws>>strw>>gkk_n_tmp>>endrw;
			stend(strw,endrw);
			rws>>strw>>gkk_n_xtmp>>endrw;
			stend(strw,endrw);
			rws>>strw>>gkk_n_ytmp>>endrw;
			stend(strw,endrw);
			rws>>strw>>gkk_n_ztmp>>endrw;
			stend(strw,endrw);
			for(int k=0;k<ngtotnod_9(i*nkpt+j);k++)
			{
				gkk(i+j+1)=gkk_n_tmp(j);
				gkk_x(i+j+1)=gkk_n_xtmp(j);
				gkk_y(i+j+1)=gkk_n_ytmp(j);
				gkk_z(i+j+1)=gkk_n_ztmp(j);
			}
			num=num+ngtotnod_9(i*nkpt+j);
		}
		ng_tot=num;
	}
	
	
	
	return true;
}
