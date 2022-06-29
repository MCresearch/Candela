#include "HB_correlation.h"
#include "gfun.h"
HB_correlation::HB_correlation(){}

HB_correlation::~HB_correlation(){}

void HB_correlation::Routine()
{
    TITLE("HB_correlation","Routine");
    ofs_running << "Calculate the time correlation function of HB." << endl;
    this->ndt = int((INPUT.msd_t)/INPUT.msd_dt);
    this->ndt1 = int((INPUT.msd_t+20)/INPUT.msd_dt);
    this->HB_start_time = new double*[INPUT.nPT];
    //this->HB_corr = new double*[INPUT.nPT];
    this->HB_nrecord = new int[INPUT.nPT]; // number of start time for each HB
    this->HB_water_index = new int*[INPUT.nPT]; // index of the HB connecting waters
    this->nHB = 0;
    this->ofs_nHB.open("nHB.txt");
    for (int iHB=0; iHB<INPUT.nPT; iHB++)
    {
        HB_start_time[iHB] = new double[this->ndt1];
        //HB_corr[iHB] = new double[ndt];
        HB_nrecord[iHB] = 0;
        HB_water_index[iHB] = new int[2];

        for (int idt=0; idt<this->ndt1; idt++)
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
        this->ofs_nHB << cel.snapshot_time << " " << this->nHB << endl;
        cel.clean();
	}//igeo
    ofs_nHB.close();
    if (INPUT.func_b == 1)
    {
        this->output();
    }
    double average_HB_lifetime = 0;
    for (int iHB=0; iHB<this->nHB; iHB++)
    {
        average_HB_lifetime += this->HB_start_time[iHB][this->HB_nrecord[iHB]-1] - this->HB_start_time[iHB][0];
    }
    average_HB_lifetime /= this->nHB;
    cout << "Average HB life time = " << average_HB_lifetime << " ps." << endl;
    ofs_running << "Average HB life time = " << average_HB_lifetime << " ps." << endl;
    return;
}

void HB_correlation::calc(CellFile &cel)
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
    for (int iwater1=0; iwater1<cel.atom[ito].na; iwater1++)
    {
        for (int idon=0; idon<water[iwater1].ndonate; idon++)
        {
            int iwater2 = water[iwater1].donateO[idon];
            int index_water_record = -1;
            for (int iHB=0; iHB<this->nHB; iHB++)
            {
                if (this->HB_water_index[iHB][0] == iwater1 and this->HB_water_index[iHB][1] == iwater2 )
                {
                    index_water_record = iHB;
                    break;
                }
            }
            if (index_water_record == -1 and this->nHB < INPUT.nPT)
            {
                this->HB_start_time[this->nHB][0] = cel.snapshot_time;
                this->HB_water_index[this->nHB][0] = iwater1;
                this->HB_water_index[this->nHB][1] = iwater2;
                this->HB_nrecord[this->nHB]++;
                //this->HB_corr[nHB][0]++;
                this->nHB++;
            }
            else if (index_water_record >= 0 and this->HB_nrecord[index_water_record] < this->ndt1)
            {
                //assert(index_water_record >= 0);
                //if (cel.snapshot_time - HB_start_time[index_water_record][0] <= INPUT.msd_t)
                //{
                this->HB_start_time[index_water_record][HB_nrecord[index_water_record]] = cel.snapshot_time;
                HB_nrecord[index_water_record]++;
                //}
            }
            else
            {
                cout << "Number of HBs reach upper boundary. Neglecting this HB." << endl;
            }
        }// end idon
    }// end iwater1
    delete[] water;
    return;
}

void HB_correlation::output()
{
    ofstream ofs_out ("HB_TCF.txt");
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