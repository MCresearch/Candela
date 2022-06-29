#include "HB_correlation2.h"
#include "gfun.h"
HB_correlation2::HB_correlation2(){}

HB_correlation2::~HB_correlation2(){}

void HB_correlation2::Routine()
{
    TITLE("HB_correlation2","Routine");
    ofs_running << "Calculate the time correlation function of HB (2nd way)." << endl;
    this->ndt = int(INPUT.msd_t/INPUT.msd_dt)+1;
    this->HB_start_time = new double**[INPUT.nPT];
    this->nslice = new int[INPUT.nPT];
    //this->HB_corr = new double**[INPUT.nPT];
    //this->HB_nrecord = new int*[INPUT.nPT]; // number of start time for each HB
    this->HB_water_index = new int*[INPUT.nPT]; // index of the HB connecting waters
    this->HB_connected = new bool[INPUT.nPT];
    this->HB_previous_connected = new bool[INPUT.nPT];
    this->nHB = 0;
    // array allocation
    for (int iHB=0; iHB<INPUT.nPT; iHB++)
    {
        this->HB_start_time[iHB] = new double*[50];
        this->nslice[iHB] = 0;
        //this->HB_corr[iHB] = new double*[20];
        //this->HB_nrecord[iHB] = new int[20];
        this->HB_water_index[iHB] = new int[2];
        this->HB_connected[iHB] = false;
        this->HB_previous_connected[iHB] = false;
        for (int islice=0; islice<50; islice++)
        {
            this->HB_start_time[iHB][islice] = new double[2];
            //this->HB_corr[iHB][islice] = new double[this->ndt];
            //this->HB_nrecord[iHB][islice] = 0;
            for(int it=0; it<2; it++)
            {
                this->HB_start_time[iHB][islice][it] = 0;
                //this->HB_corr[iHB][islice][it] = 0;
            }
        }
        for (int ia=0; ia<2; ia++)
        {
            HB_water_index[iHB][ia] = -1;
        }
    }
    // end array allocation
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
	}//igeo

    this->output();
    return;
}

void HB_correlation2::calc(CellFile &cel)
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
            if (index_water_record == -1)
            {
                this->HB_water_index[this->nHB][0] = iwater1;
                this->HB_water_index[this->nHB][1] = iwater2;
                this->HB_start_time[this->nHB][0][0] = cel.snapshot_time;
                //this->HB_nrecord[this->nHB][0]++;
                //this->HB_corr[this->nHB][0][0]++;
                this->HB_connected[this->nHB] = true;
                this->nHB++;
            }
            if (index_water_record >= 0)
            {
                if (!this->HB_connected[index_water_record] and this->nslice[index_water_record]<20)
                {
                    this->HB_start_time[index_water_record][this->nslice[index_water_record]][0] = cel.snapshot_time;
                    this->HB_connected[index_water_record] = true;
                    //this->HB_corr[index_water_record][this->nslice[index_water_record]][0]++;
                    //this->HB_nrecord[index_water_record][this->nslice[index_water_record]]++;
                }
                //else if (this->HB_connected[index_water_record])
                // {
                //    if (cel.snapshot_time - this->HB_start_time[index_water_record][this->nslice[index_water_record]][0] < INPUT.msd_t)
                //    {
                //        this->HB_connected[index_water_record] = true;
                //    }
                //}
            }
        }
    }
    for (int iHB=0; iHB<this->nHB; iHB++)
    {
        if (this->HB_connected[iHB])
        {
            int water_index1 = this->HB_water_index[iHB][0];
            int water_index2 = this->HB_water_index[iHB][1];
            bool HB_con = false;
            for (int idon=0; idon<water[water_index1].ndonate; idon++)
            {
                if(water[water_index1].donateO[idon] == water_index2)
                {
                    HB_con = true;
                    break;
                }
            }
            if (!HB_con)
            {
                //cout << "iHB = " << iHB << " not connected." << endl;
                this->HB_connected[iHB] = false;
                this->HB_start_time[iHB][this->nslice[iHB]][1] = cel.snapshot_time-INPUT.msd_dt;
                this->nslice[iHB]++;
                //cout << this->nslice[iHB]++ << endl;
            }
        }
        //cout << "iHB = " << iHB << ", nHB = " << this->nHB << endl;
    }
}

void HB_correlation2::output()
{
    double* TCF = new double[this->ndt];
    ofstream ofs("HB_TCF.txt");
    for (int idt=0; idt<this->ndt; idt++)
    {
        TCF[idt] = 0;
    }
    for (int iHB=0; iHB<this->nHB; iHB++)
    {
        for (int islice=0; islice<this->nslice[iHB]; islice++)
        {
            double upper_index = round((this->HB_start_time[iHB][islice][1] - this->HB_start_time[iHB][islice][0])/INPUT.msd_dt);
            double portion = 1/upper_index;
            int upper_index1 = int(upper_index);
            for (int it=0; it<upper_index1; it++)
            {
                TCF[it]+=portion*(upper_index-it);
            }
        }
    }
    double total_nslice = TCF[0];
    for (int idt=0; idt<this->ndt; idt++)
    {
        TCF[idt] /= total_nslice;
        ofs << idt*INPUT.msd_dt << " " << TCF[idt] << endl;
    }
    ofs.close();
    return;
}