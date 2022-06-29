#include "HB_stat4.h"
#include "dist2.h"
#include "gfun.h"
HB_stat4::HB_stat4(){}
HB_stat4::~HB_stat4(){}

void HB_stat4::Routine()
{
    TITLE("HB_stat4","Routine");
    ofs_running << "HB orientation statistics." << endl;
    int count_geometry_number = 0;
    this->ndouble_donate = 0;
    this->nHB_jump = 0;
    this->ofs_DD.open("double_donor.txt");
    this->nspace = 20;
    this->ofs_DD << setw(nspace) << "start_time" << setw(nspace) << "end_time"
    << setw(nspace) << "central_O_index" << setw(nspace) << "donated_H_index" << 
     setw(nspace) << "Oa_index" << setw(nspace) << "Ob_index" <<  setw(nspace) << 
     "previous_donated_O" << setw(nspace) << "following_donate_O" << setw(nspace) << 
     "HBjump"  << setw(nspace) << "Hconnected" << endl;
    this->ofs = new ofstream[INPUT.natom_new];
    this->previous_donateO = new int[INPUT.natom2];
    this->previous_doubledonateO = new int*[INPUT.natom2];
    this->now_donateO = new int[INPUT.natom2];
    this->previous_double_donate = new bool[INPUT.natom2];
    this->double_donate = new bool[INPUT.natom2];
    this->DD_start_time = new double[INPUT.natom2];
    this->DD_HB = new string[INPUT.natom2];
    for (int ia=0; ia<INPUT.natom_new; ia++)
    {
        this->ofs[ia].open("water_donate_" + to_string(ia) + ".txt");
    }
    for (int ia=0; ia<INPUT.natom2; ia++)
    {
        this->previous_donateO[ia] = -1;
        this->now_donateO[ia] = -1;
        this->previous_double_donate[ia] = false;
        this->double_donate[ia] = false;
        this->DD_start_time[ia] = -1.0;
        this->previous_doubledonateO[ia] = new int[2];
        this->previous_doubledonateO[ia][0] = -1;
        this->previous_doubledonateO[ia][1] = -1;
        this->DD_HB[ia] = "None";
    }
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
    for (int ia=0; ia<INPUT.natom_new; ia++)
    {
        ofs[ia].close();
    }
    cout << "Total number of double donors = " << this->ndouble_donate << endl;
    cout << "Number of HB jumps in double donors = " << this->nHB_jump << endl;
    ofs_running << "Total number of double donors = " << this->ndouble_donate << endl;
    ofs_running << "Number of HB jumps in double donors = " << this->nHB_jump << endl;
    this->ofs_DD.close();
    return;
}

void HB_stat4::calc(CellFile &cel)
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
    int** donateO = new int*[2];
    for (int i = 0; i<2; i++)
    {
        donateO[i] = new int[5];
        for (int j = 0; j < 5; j++)
        {
            donateO[i][j] = -1;
        }
    }
    int* ndonate = new int[2];
    for (int iwater=0; iwater<cel.atom[ito].na; iwater++)
    {
        for (int i = 0; i<2; i++)
        {
            ndonate[i] = 0;
            for (int j = 0; j < 5; j++)
            {
                donateO[i][j] = -1;
            }
        }
        for (int idonate=0; idonate<water[iwater].ndonate; idonate++)
        {
            for (int iH=0; iH<2; iH++)
            {
                if ( water[iwater].donateH[idonate] == water[iwater].indexH[iH] )
                {
                    donateO[iH][ndonate[iH]] = water[iwater].donateO[idonate];
                    ndonate[iH]++;
                }
            }
        }
        for (int iH=0; iH < 2; iH++)
        {
            int Hindex = water[iwater].indexH[iH];
            if (ndonate[iH]>=2)
            {
                this->double_donate[Hindex] = true;
                this->previous_doubledonateO[Hindex][0] = donateO[iH][0];
                this->previous_doubledonateO[Hindex][1] = donateO[iH][1];
            }
            else
            {
                this->double_donate[Hindex] = false;
            }
            if (!this->double_donate[Hindex])
            {
                this->now_donateO[Hindex] = donateO[iH][0];
            }
            if (this->previous_double_donate[Hindex] == false and 
            this->double_donate[Hindex] == true)
            {
                this->DD_start_time[Hindex] = cel.snapshot_time;
                bool DD_HB_flag = false;
                for (int iacc=0; iacc<water[donateO[iH][0]].naccept; iacc++)
                {
                    if (water[donateO[iH][0]].acceptO[iacc] == donateO[iH][1]) DD_HB_flag = true;
                    
                }
                for (int idon=0; idon<water[donateO[iH][0]].ndonate; idon++)
                {
                    if (water[donateO[iH][0]].acceptO[idon] == donateO[iH][1]) DD_HB_flag = true;
                }
                if(DD_HB_flag) this->DD_HB[Hindex] = "Connected";
                if(!DD_HB_flag) this->DD_HB[Hindex] = "Disconnected";
            }

            if (this->previous_double_donate[Hindex] == true and 
            this->double_donate[Hindex] == false)
            {
                this->ndouble_donate++;
                string jump = "False";
                if (this->now_donateO[Hindex] != this->previous_donateO[Hindex] and
                this->previous_donateO[Hindex] != -1)
                {
                    this->nHB_jump++;
                    jump = "True";
                }
                this->ofs_DD << setw(nspace) << DD_start_time[Hindex] << setw(nspace) << cel.snapshot_time << 
                setw(nspace) << iwater << setw(nspace) << Hindex << setw(nspace) <<  this->previous_doubledonateO[Hindex][0] << setw(nspace) << 
                this->previous_doubledonateO[Hindex][1] << setw(nspace) << this->previous_donateO[Hindex] << setw(nspace) << 
                this->now_donateO[Hindex] << setw(nspace) << jump << setw(nspace) << this->DD_HB[Hindex] << endl;
            }
            this->previous_donateO[Hindex] = this->now_donateO[Hindex];
            this->previous_double_donate[Hindex] = this->double_donate[Hindex];

            if(iwater < INPUT.natom_new)
            {
                this->ofs[iwater] << setprecision(10) << setw(20) << cel.snapshot_index << setw(20) << cel.snapshot_time << setw(20) << iwater
                << setw(20) << Hindex << setw(20) << ndonate[iH];
                //cout << setprecision(10) << setw(20) << cel.snapshot_index << setw(20) << cel.snapshot_time << setw(20) << iwater
                // << setw(20) << water[iwater].indexH[iH] << setw(20) << ndonate[iH] << endl;
                for (int isdon = 0; isdon<5; isdon++)
                {
                    this->ofs[iwater] << setw(20) << donateO[iH][isdon];
                }
                this->ofs[iwater] << endl;
            }
        }
    }
    for (int i = 0; i<2; i++)
    {
        delete[] donateO[i];
    }
    delete[] donateO;
    delete[] ndonate;
    delete[] water;
    return;
}