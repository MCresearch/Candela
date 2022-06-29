#include "HB_break.h"

HB_break::HB_break(){}

HB_break::~HB_break(){}

void HB_break::Routine()
{
    TITLE("HB_break","Routine");
    ofs_running << "Count the reasons for HB break-up." << endl;

    this->last_donate_index = new int*[INPUT.natom1];
    this->last_ndonate = new int[INPUT.natom1];
    this->last_donate_H = new int*[INPUT.natom1];
    for (int ia=0; ia<INPUT.natom1; ia++)
    {
        this->last_donate_index[ia] = new int[5];
        this->last_donate_H[ia] = new int[5];
        this->last_ndonate[ia] = 0;
        for (int idon=0; idon<5; idon++)
        {
            this->last_donate_index[ia][idon] = -1;
            this->last_donate_H[ia][idon] = -1;
        }
    } 
    this->break_reason_count = new double[3];
    this->form_reason_count = new double[3];
    for (int i=0; i<3; i++)
    {
        this->break_reason_count[i] = 0;
        this->form_reason_count[i] = 0;
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
	}//igeo    
    int sum_reason = 0;
    sum_reason = this->break_reason_count[0]+this->break_reason_count[1]+this->break_reason_count[2];
    ofs_running << "Rcut exceeded: " << this->break_reason_count[0]/sum_reason << endl;
    ofs_running << "Acut exceeded: " << this->break_reason_count[1]/sum_reason << endl;
    ofs_running << "Both Rcut and Acut exceeded: " << this->break_reason_count[2]/sum_reason << endl;
    return;
}


void HB_break::calc(CellFile &cel)
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

    for (int iwater=0; iwater<cel.atom[ito].na; iwater++)
    {
        if (this->last_ndonate[iwater] > 0)
        {
            for (int idon=0; idon<this->last_ndonate[iwater]; idon++)
            {
                if (!in_list(this->last_donate_index[iwater][idon], water[iwater].donateO, water[iwater].ndonate))
                {
                    double dis_oo = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[water[iwater].donateO[idon]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                    double ang_ooh = HBs::angle(cel, cel.atom[ito].pos[this->last_donate_index[iwater][idon]], cel.atom[ito].pos[iwater], cel.atom[ith].pos[this->last_donate_H[iwater][idon]]);
                    if (dis_oo > INPUT.rcut_oo and ang_ooh > INPUT.acut_hoo)
                    {
                        this->break_reason_count[2]++;
                    }
                    else if (dis_oo > INPUT.rcut_oo)
                    {
                        this->break_reason_count[0]++;
                    }
                    else if (ang_ooh > INPUT.acut_hoo)
                    {
                        this->break_reason_count[1]++;
                    }
                }
            }
        }
        for (int idon=0; idon<this->last_ndonate[iwater]; idon++) 
        {
            this->last_donate_H[iwater][idon] = -1;
            this->last_donate_index[iwater][idon] = -1;
        }
        this->last_ndonate[iwater] = water[iwater].ndonate;
        for (int idon=0; idon<this->last_ndonate[iwater]; idon++)
        {
            this->last_donate_H[iwater][idon] = water[iwater].donateH[idon];
            this->last_donate_index[iwater][idon] = water[iwater].donateO[idon];
        }

    }
    delete[] water;
}

bool HB_break::in_list(int &target, int* list, int &list_length)
{
    for (int index=0; index<list_length; index++)
    {
        if (list[index] == target)
        {
            return true;
        }
    }
    return false;
}