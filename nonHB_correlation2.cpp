#include "nonHB_correlation2.h"

nonHB_correlation2::nonHB_correlation2(){}

nonHB_correlation2::~nonHB_correlation2(){}

void nonHB_correlation2::Routine()
{
    TITLE("HB_correlation","Routine");
    ofs_running << "Calculate the time correlation function of HB." << endl;
    this->ndt = int(INPUT.msd_t/INPUT.msd_dt)+1;
    this->ndt1 = int((INPUT.msd_t+20)/INPUT.msd_dt);
    this->HB_start_time = new double*[INPUT.nPT];
    //this->HB_corr = new double*[INPUT.nPT];
    this->HB_nrecord = new int[INPUT.nPT]; // number of start time for each HB
    this->HB_water_index = new int*[INPUT.nPT]; // index of the HB connecting waters
    this->nHB = 0;
    this->ofs_nHB.open("nHB.txt");
    for (int iHB=0; iHB<INPUT.nPT; iHB++)
    {
        HB_start_time[iHB] = new double[ndt1];
        //HB_corr[iHB] = new double[ndt];
        HB_nrecord[iHB] = 0;
        HB_water_index[iHB] = new int[2];

        for (int idt=0; idt<ndt1; idt++)
        {
            HB_start_time[iHB][idt] = -1;
        }
        for (int iwater=0; iwater<2; iwater++)
        {
            HB_water_index[iHB][iwater] = -1;
        }
    }

    int count_geometry_number = 0;
    for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		// cel : input geometry file
		CellFile cel;

		//cout << " igeo=" << igeo << " igeo%INPUT.geo_interval=" << igeo%INPUT.geo_interval << endl;
        if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) 
		{
			cel.clean(); // renxi added 20200614
			continue;
		}
		++count_geometry_number;
		cout << "snapshot " << igeo << endl;
        this->calc(cel);
        cel.clean();
        this->ofs_nHB << cel.snapshot_time << " " << this->nHB << endl;
	}//igeo
    if (INPUT.func_b == 1)
    {
        this->output();
    }
    this->ofs_nHB.close();
    return;
}

void nonHB_correlation2::calc(CellFile &cel)
{
    int ito=-1;
	int ith=-1;
	int itc=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="C") itc=it;
	}
	if(INPUT.ntype==2){ assert(ito>=0); assert(ith>=0);}
    
    Water *water = new Water[cel.atom[ito].na];
    Water::nions = 0;
    HBs::setup_water(cel, water);
    for (int iwater1=0; iwater1<cel.atom[ito].na-1; iwater1++)
    {
        for (int iwater2=iwater1+1; iwater2<cel.atom[ito].na; iwater2++)
        {
            double dist_oo = distance(cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
            if (dist_oo < INPUT.rcut_oo)
            {
                int index_of_HB = -1;
                for (int iHB=0; iHB<this->nHB; iHB++)
                {
                    if (this->HB_water_index[iHB][0] == iwater1 and this->HB_water_index[iHB][1] == iwater2)
                    {
                        index_of_HB = iHB;
                        break;
                    }
                }
                if (index_of_HB >= 0)
                {
                    if (this->HB_nrecord[index_of_HB] < this->ndt1)
                    {
                        this->HB_start_time[index_of_HB][this->HB_nrecord[index_of_HB]] = cel.snapshot_time;
                        this->HB_nrecord[index_of_HB]++;
                    }
                }
                else if (this->nHB < INPUT.nPT)
                {
                    this->HB_water_index[this->nHB][0] = iwater1;
                    this->HB_water_index[this->nHB][1] = iwater2;
                    this->HB_start_time[this->nHB][0] = cel.snapshot_time;
                    this->HB_nrecord[this->nHB] = 1;
                    this->nHB++;
                }
            }
        }
    }
    delete[] water;
    return;
}

void nonHB_correlation2::output()
{
    ofstream ofs_out ("nonHB_TCF.txt");
    double* tcf = new double[this->ndt];
    double* tcf_total = new double[this->ndt];
    for (int idt=0; idt<this->ndt; idt++) tcf_total[idt] = 0;
    for (int iHB=0; iHB<this->nHB; iHB++)
    {
        for (int idt=0; idt<this->ndt; idt++) tcf[idt] = 0;
        for (int irecord1=0; irecord1<this->HB_nrecord[iHB]; irecord1++)
        {
            for (int irecord2=irecord1; irecord2<this->HB_nrecord[iHB]; irecord2++)
            {
                int which_time = int(round((this->HB_start_time[iHB][irecord2] - this->HB_start_time[iHB][irecord1])/INPUT.msd_dt));
                assert(which_time >= 0);
                if (which_time < this->ndt)
                {
                    tcf[which_time]++;
                }
            }
        }
        for (int idt=0; idt<this->ndt; idt++) 
        {
            tcf[idt] /= this->HB_nrecord[iHB];
            tcf_total[idt] += tcf[idt];
        }
    }
    for (int idt=0; idt<this->ndt; idt++)
    {
        tcf_total[idt] /= this->nHB;
        ofs_out << idt*INPUT.msd_dt << " " << tcf_total[idt] << endl;
    }
    ofs_out.close();
    return;
}