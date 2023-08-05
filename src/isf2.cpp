#include "isf2.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include "cellFile.h"
#include "input.h"
#include "vec3.h"
void ISF2::Routine()
{
	TITLE("ISF2","Routine");
	cal();
	return;
}
void ISF2::cal()
{
		int allcount=-1;
		int count=0;
        for(int m1=0;m1<INPUT.isf_ngx;m1++)
        {
          for(int m2=0;m2<INPUT.isf_ngy;m2++)
          {
            for(int m3=0;m3<INPUT.isf_ngz;m3++)
            {
        		if(!targetk(m1,m2,m3))
        		{
            		continue;
        		}
				allcount++;
#ifdef __MPI
				if((allcount)%NPROC==RANK)
#endif
				count++;
			}
		  }
		}
		assert(count<1000);
		if(count == 0)
		{
			cout<<"Find no q! Please change target_q!"<<endl;
			exit(0);
		}
		Vector3<int>* mg=new Vector3<int>[count];
		int *pre=new int[count];
		allcount=count=-1;
        for(int m1=0;m1<INPUT.isf_ngx;m1++)
        {
          for(int m2=0;m2<INPUT.isf_ngy;m2++)
          {
            for(int m3=0;m3<INPUT.isf_ngz;m3++)
            {
                if(!targetk(m1,m2,m3))
                {
                        continue;
                }
				allcount++;
#ifdef __MPI
				if((allcount)%NPROC==RANK)
#endif
				{
				count++;
	        	if((m1==0&&m2==0)||(m1==0&&m3==0)||(m2==0&&m3==0)) pre[count]=1;
                else if(m1==0||m2==0||m3==0) pre[count]=2;
                else pre[count]=4;
        		cout<<'('<<m1<<','<<m2<<','<<m3<<") is used."<<endl;
				mg[count].x=m1;
				mg[count].y=m2;
				mg[count].z=m3;
				}
			}
		  }
		}
		count++;
		if(RANK==0)		cout<<"Total: "<<allcount+1<<endl;
		double *rhok=new double[count*INPUT.isf_nt1*8];
        double *ddcf=new double[count*INPUT.isf_nt1];
        int *index=new int[count*INPUT.isf_nt1];	
		for(int i=0;i<count*INPUT.isf_nt1;i++)
		{
          	ddcf[i]=0;
           	index[i]=0;
		}
		//ignore geos before geo_1
		geoignore();
    	//get ddcf(density-density correlation function)
		for(int tt=0;tt<INPUT.isf_nt1;tt++)
		{
			if(!getrhok(tt,mg,tt,rhok,count))
			{
				cout<<"Error in loading geo!"<<endl;
			}
		}
		
		for(int tt=0;tt<INPUT.isf_nt1;tt++)
		{
			int iT=0;
			while(iT+tt<INPUT.isf_nt1)
            {
                for(int mm=0;mm<count;mm++)
				{
					for(int i=0;i<8;i++)
                	{
                		ddcf[tt*count+mm]+=pre[mm]*rhok[(tt+iT)*count*8+mm*8+i]*rhok[iT*count*8+mm*8+i];
                	}
					index[tt*count+mm]+=pre[mm];
				}
                iT++;
          	}
		}
		for(int TT=0;TT<INPUT.isf_nt2;TT++)
		{
			if(getigeo(TT+INPUT.isf_nt1)>INPUT.geo_2) continue;
        	if(!getrhok(TT+INPUT.isf_nt1,mg,TT%INPUT.isf_nt1,rhok,count))
        	{
            	cout<<"erro in loading geo"<<endl;
            	return;
        	}
            for(int tt=0;tt<INPUT.isf_nt1;tt++)
            {
            	for(int mm=0;mm<count;mm++)
				{
					for(int i=0;i<8;i++)
        			{     
						ddcf[tt*count+mm]+=pre[mm]*rhok[(TT%INPUT.isf_nt1)*count*8+mm*8+i]*rhok[((TT+INPUT.isf_nt1-tt)%INPUT.isf_nt1)*count*8+mm*8+i];
                	}
                	index[tt*count+mm]+=pre[mm];
				}
            }
        }
		delete[] mg;
		delete[] pre;
//converge different k
double* cv_ddcf=new double[INPUT.isf_nt1];
int* cv_index=new int[INPUT.isf_nt1];
ZEROS(cv_ddcf,INPUT.isf_nt1);
ZEROS(cv_index,INPUT.isf_nt1);
for(int tt=0;tt<INPUT.isf_nt1;tt++)
{
	for(int mm=0;mm<count;mm++)
	{
		cv_ddcf[tt]+=ddcf[tt*count+mm];
		cv_index[tt]+=index[tt*count+mm];
	}
}

#ifdef __MPI
	double* f_ddcf=new double[INPUT.isf_nt1];
    int* f_index=new int[INPUT.isf_nt1];
    MPI_Reduce(cv_ddcf,f_ddcf,INPUT.isf_nt1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(cv_index,f_index,INPUT.isf_nt1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	writeisf(f_ddcf,f_index);
	delete []f_ddcf;
	delete []f_index;
#else
	writeisf(cv_ddcf,cv_index);
#endif
	//clear
	delete []rhok;
    delete []ddcf;
    delete []cv_ddcf;
    delete []index;
    delete []cv_index;
	return;
}

bool targetk(int m1,int m2,int m3)
{
        //borhr=0.529177A
		double mol=sqrt(pow((m1*INPUT.isf_dgx),2)+pow((m2*INPUT.isf_dgy),2)+pow((m3*INPUT.isf_dgz),2));
        if(abs(mol-INPUT.isf_target_q)<0.005) return true;
        else return false;
}
double gettime(int tt)
{
        return INPUT.geo_interval*INPUT.dt_snapshots*tt;
}
int getigeo(int tt)
{
	return INPUT.geo_interval*tt+INPUT.geo_1;
}
template<class T>
bool getrhok(int ntt,Vector3<T>* mg,int tt,double *rho,int &count)
{
	CellFile cel_in;

	cout<<"igeo="<<getigeo(ntt)<<endl;

	stringstream ss; ss << getigeo(ntt);
        cel_in.file_name = ss.str();
	//read in geometry
	if(!CellFile::ReadGeometry(cel_in)) return false;
	
	
        Vector3<double> r;
        Vector3<double> cosp,sinp;
        Vector3<double>  phi;
        double* tem=new double[8];
	for(int mm=0;mm<count;mm++)
	{
		for(int i=0;i<8;i++) tem[i]=0;
        for(int ia=0;ia<INPUT.natom;ia++)
        {
			r=cel_in.atom[0].pos[ia];
			phi.x=mg[mm].x*INPUT.isf_dgx*r.x;
            phi.y=mg[mm].y*INPUT.isf_dgy*r.y;
            phi.z=mg[mm].z*INPUT.isf_dgz*r.z;
            cosp=cos(phi);
            sinp=sin(phi);
            tem[0]+=cosp.x*cosp.y*cosp.z;
            tem[1]+=cosp.x*cosp.y*sinp.z;
            tem[2]+=cosp.x*sinp.y*cosp.z;
            tem[3]+=cosp.x*sinp.y*sinp.z;
            tem[4]+=sinp.x*cosp.y*cosp.z;
            tem[5]+=sinp.x*cosp.y*sinp.z;
            tem[6]+=sinp.x*sinp.y*cosp.z;
            tem[7]+=sinp.x*sinp.y*sinp.z;
        }
        for (int i=0;i<8;i++)
        {
                rho[tt*count*8+mm*8+i]=tem[i];
        }
	}
        delete []tem;
	cel_in.clean();
        return true;
}
void geoignore()
{
	for(int i=0;i<INPUT.geo_1;i+=INPUT.geo_interval)
	{
		if (INPUT.geo_in_type=="PROFESS") break;
		cout<<"ignore_igeo: "<<i<<endl;
		CellFile cel;
		CellFile::ReadGeometry(cel);
		cel.clean();
	}
	return;
}
void writeisf(double *f_ddcf,int *f_index)
{
    if(RANK==0)
	{
			for(int tt=0;tt<INPUT.isf_nt1;tt++)
        	{
           	     	f_ddcf[tt]=f_ddcf[tt]/INPUT.natom/f_index[tt];
        	}
        	//print results
        	ofstream ofs(INPUT.isf_outfile.c_str());
        	ofs<<"#time ISF"<<endl;
        	for (int tt=0;tt<INPUT.isf_nt1;tt++)
        	{
                	ofs<<gettime(tt)<<' '<<f_ddcf[tt]<<endl;
        	}
        	ofs.close();
	}
	return;
}
