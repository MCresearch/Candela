#include "nonHB_correlation.h"
#include "gfun.h"
nonHB_correlation::nonHB_correlation(){}

nonHB_correlation::~nonHB_correlation(){}

void nonHB_correlation::Routine()
{
    TITLE("nonHB_correlation","Routine");
    ofs_running << "Calculate the time correlation function of non-HB inner shell TCF ." << endl;
    this->ndt = int(INPUT.msd_t/INPUT.msd_dt)+1;
    this->nonHB_time = new double**[INPUT.nPT];
    this->nslice = new int[INPUT.nPT];
    //this->HB_corr = new double**[INPUT.nPT];
    //this->HB_nrecord = new int*[INPUT.nPT]; // number of start time for each HB
    this->HB_water_index = new int*[INPUT.nPT]; // index of the HB connecting waters
    this->inshell = new bool[INPUT.nPT];
    this->previous_inshell = new bool[INPUT.nPT];
    this->nHB = 0;
    // array allocation
    for (int iHB=0; iHB<INPUT.nPT; iHB++)
    {
        this->nonHB_time[iHB] = new double*[200];
        this->nslice[iHB] = 0;
        //this->HB_corr[iHB] = new double*[20];
        //this->HB_nrecord[iHB] = new int[20];
        this->HB_water_index[iHB] = new int[2];
        this->inshell[iHB] = false;
        this->previous_inshell[iHB] = false;
        for (int islice=0; islice<200; islice++)
        {
            this->nonHB_time[iHB][islice] = new double[2];
            //this->HB_corr[iHB][islice] = new double[this->ndt];
            //this->HB_nrecord[iHB][islice] = 0;
            for(int it=0; it<2; it++)
            {
                this->nonHB_time[iHB][islice][it] = 0;
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
}


void nonHB_correlation::calc(CellFile &cel)
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
            bool judge = ((INPUT.func == 1 and dist_oo <= INPUT.rcut_oo) or (INPUT.func == 2 and dist_oo <= INPUT.rcut_oo and !this->HB_bonded(water, iwater1, iwater2)));
            if (judge)
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
                if (index_of_HB < 0)
                {
                    this->HB_water_index[this->nHB][0] = iwater1;
                    this->HB_water_index[this->nHB][1] = iwater2;
                    this->nonHB_time[this->nHB][0][0] = cel.snapshot_time;
                    //this->nslice[this->nHB]++;
                    this->inshell[this->nHB] = true;
                    this->previous_inshell[this->nHB] = true;
                    this->nHB++;
                }
                else
                {
                    if (!this->previous_inshell[index_of_HB])
                    {
                        this->nonHB_time[index_of_HB][this->nslice[index_of_HB]][0] = cel.snapshot_time;
                        this->inshell[index_of_HB] = true;
                        //this->nslice[index_of_HB]++;
                        this->previous_inshell[index_of_HB] = true;
                    }
                }
            }
        }
    }
    for (int iHB=0; iHB<this->nHB; iHB++)
    {
        int iwater1=this->HB_water_index[iHB][0];
        int iwater2=this->HB_water_index[iHB][1];
        double dist_oo = distance(cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
        if (INPUT.func == 1)
        {
            bool judge = (dist_oo <= INPUT.rcut_oo);
            if (!judge)
            {
                if (this->previous_inshell[iHB])
                {
                    this->nonHB_time[iHB][this->nslice[iHB]][1] = cel.snapshot_time;
                    this->nslice[iHB]++;
                    this->previous_inshell[iHB] = false;
                    this->inshell[iHB] = false;
                }
            }
        }
        else if (INPUT.func == 2)
        {
            bool judge = (this->HB_bonded(water, iwater1, iwater2));
            if (judge)
            {
                if (this->previous_inshell[iHB])
                {
                    this->nonHB_time[iHB][this->nslice[iHB]][1] = cel.snapshot_time;
                    this->nslice[iHB]++;
                    this->previous_inshell[iHB] = false;
                    this->inshell[iHB] = false;
                }
            }
        }
    }
    delete[] water;
    return;
}

bool nonHB_correlation::HB_bonded(Water* water, int &iwater1, int &iwater2)
{
    for (int idon=0; idon<water[iwater1].ndonate; idon++)
    {
        //cout << idon << " " << water[iwater1].donateO[idon] << endl;
        if (water[iwater1].donateO[idon] == iwater2)
        {
            return true;
        }
    }
    for (int iacc=0; iacc<water[iwater1].naccept; iacc++)
    {
        //cout << iacc << " " << water[iwater1].acceptO[iacc] << endl;
        if (water[iwater1].acceptO[iacc] == iwater2)
        {
            return true;
        }
    }
    return false;
}

void nonHB_correlation::output()
{
    double* TCF = new double[this->ndt];
    ofstream ofs("nonHB_TCF.txt");
    for (int idt=0; idt<this->ndt; idt++)
    {
        TCF[idt] = 0;
    }
    for (int iHB=0; iHB<this->nHB; iHB++)
    {
        for (int islice=0; islice<this->nslice[iHB]; islice++)
        {
            double upper_index = round((this->nonHB_time[iHB][islice][1] - this->nonHB_time[iHB][islice][0])/INPUT.msd_dt);
            double portion = 1/upper_index;
            int upper_index1 = int(upper_index);
            upper_index1 = min(upper_index1, this->ndt-1);
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