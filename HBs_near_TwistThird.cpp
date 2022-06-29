#include "HBs_near_TwistThird.h"
#include "gfun.h"
HBs_near_TwistThird::HBs_near_TwistThird(){}

HBs_near_TwistThird::~HBs_near_TwistThird(){}

void HBs_near_TwistThird::Routine()
{
    TITLE("HBs_near_TwistThird","Routine");
    ofs_running << "Check whether the third neighbor is HBed to the central O atom." << endl;

    int count_geometry_number = 0;
    this->water_index = new int*[INPUT.nPT];
    this-> start_end_time = new double*[INPUT.nPT];
    this->HB_before_amid_after = new bool*[INPUT.nPT];
    for (int iTT=0; iTT<INPUT.nPT; iTT++)
    {
        water_index[iTT] = new int[2];
        start_end_time[iTT] = new double[2];
        HB_before_amid_after[iTT] = new bool[3];
        for (int i=0; i<3; i++) HB_before_amid_after[iTT][i] = false;
    }
    this->allocate_index_time();
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
    int nHB_before = 0;
    int nHB_amid = 0;
    int nHB_after = 0;
    for (int iTT=0; iTT<INPUT.nPT; iTT++)
    {
        if (this->HB_before_amid_after[iTT][0]) nHB_before++;
        if (this->HB_before_amid_after[iTT][1]) nHB_amid++;
        if (this->HB_before_amid_after[iTT][2]) nHB_after++;
    }
    cout << "Percentage of HBed twisted third neighbor in " << INPUT.msd_dt0 << " ps before: " << nHB_before << "/" << INPUT.nPT << endl;
    cout << "Percentage of mixed HBed and twisted third neighbor " << nHB_amid << "/" << INPUT.nPT << endl;
    cout << "Percentage of HBed twisted third neighbor in " << INPUT.msd_dt0 << " ps after: " << nHB_after << "/" << INPUT.nPT << endl;
    ofstream ofs_before_amid_after("HB_before_amid_after.txt");
    for (int iTT=0; iTT<INPUT.nPT; iTT++)
    {
        for (int i=0; i<3; i++)
        {
            ofs_before_amid_after << this->HB_before_amid_after[iTT][i] << " ";
        }
        ofs_before_amid_after << endl;
    }
    ofs_before_amid_after.close();
    return;
}

void HBs_near_TwistThird::allocate_index_time()
{
    ifstream ifs_record("twist_third_neighbor_record.txt");
    string useless;
    for (int iTT=0; iTT<INPUT.nPT; iTT++)
    {
        ifs_record >> this->start_end_time[iTT][0] >> this->start_end_time[iTT][1] >> this->water_index[iTT][0] >> useless >> useless >> this->water_index[iTT][1];
    }
    return;
}

void HBs_near_TwistThird::calc(CellFile &cel)
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
    for (int iTT=0; iTT<INPUT.nPT; iTT++)
    {
        if (cel.snapshot_time - this->start_end_time[iTT][0] < INPUT.msd_dt0)
        {
            if (Hbonded(water, this->water_index[iTT][0], this->water_index[iTT][1]))
            {
                this->HB_before_amid_after[iTT][0] = true;
            }
        }
        if (cel.snapshot_time >= this->start_end_time[iTT][0] and cel.snapshot_time <= this->start_end_time[iTT][1])
        {
            if (Hbonded(water, this->water_index[iTT][0], this->water_index[iTT][1]))
            {
                this->HB_before_amid_after[iTT][1] = true;
            }
        }
        if (this->start_end_time[iTT][1] - cel.snapshot_time < INPUT.msd_dt0)
        {
            if (Hbonded(water, this->water_index[iTT][0], this->water_index[iTT][1]))
            {
                this->HB_before_amid_after[iTT][2] = true;
            }
        }
    }
    delete[] water;
    return;
}

bool HBs_near_TwistThird::Hbonded(Water* water, int &iwater1, int &iwater2)
{
    bool HB = false;
    for (int idon=0; idon < water[iwater1].ndonate; idon++)
    {
        if (water[iwater1].donateO[idon] == iwater2)
        {
            HB = true;
            return HB;
        }
    }
    for (int iacc=0; iacc < water[iwater1].naccept; iacc++)
    {
        if(water[iwater1].acceptO[iacc] == iwater2)
        {
            HB = true;
            return HB;
        }
    }
    return HB;
}