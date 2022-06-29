#include "HBs_near_DoubleDonor.h"
#include "gfun.h"
#include "HBs_near_AngularJump.h"
HBs_near_DoubleDonor::HBs_near_DoubleDonor(){}

HBs_near_DoubleDonor::~HBs_near_DoubleDonor(){}

void HBs_near_DoubleDonor::Routine()
{
    TITLE("HB_stat3","Routine");
    ofs_running << "HB reorientation chrononical statistics." << endl;
    
    assert(INPUT.nPT > 0); // nPT is used to count number of angular jumps.
    assert(INPUT.stay_tmax > 0); // time before and after Angular jump counted.
    assert(INPUT.stay_dt > 0); 
    assert(INPUT.relax_time > 0); // the time before AJ happens to document Oa and Ob index.

    this->ndt = int(2*INPUT.stay_tmax/INPUT.stay_dt)+5;
    //int ndr = int(INPUT.rcut/INPUT.dr);
    cout << "ndt = " << this->ndt << endl;
    DD = new DoubleDonor[INPUT.nPT];
    DoubleDonor::setup_DD(INPUT.nPT, this->DD);

    this->OOO_angle = new double*[INPUT.nPT];
    this->OH_central_angle = new double*[INPUT.nPT];
    this->OH_plane_angle = new double*[INPUT.nPT];
    this->OaOH_angle = new double*[INPUT.nPT];
    this->ObOH_angle = new double*[INPUT.nPT];
    this->OaO_distance = new double*[INPUT.nPT];
    this->ObO_distance = new double*[INPUT.nPT];

    this->central_anthH_donate = new double*[INPUT.nPT];
    this->central_accept = new double*[INPUT.nPT];
    this->donOa_donate = new double*[INPUT.nPT];
    this->donOb_donate = new double*[INPUT.nPT];
    this->donOa_accept = new double*[INPUT.nPT];
    this->donOb_accept = new double*[INPUT.nPT];
    for (int iAJ = 0; iAJ < INPUT.nPT; iAJ++)
    {
        
        this->OH_central_angle[iAJ] = new double[ndt];
        this->OH_plane_angle[iAJ] = new double[ndt];
        this->OaOH_angle[iAJ] = new double[ndt];
        this->ObOH_angle[iAJ] = new double[ndt];
        this->OaO_distance[iAJ] = new double[ndt];
        this->ObO_distance[iAJ] = new double[ndt];
        this->OOO_angle[iAJ] = new double[ndt];
        this->central_anthH_donate[iAJ] = new double[ndt];
        this->central_accept[iAJ] = new double[ndt];
        this->donOa_donate[iAJ] = new double[ndt];
        this->donOb_donate[iAJ] = new double[ndt];
        this->donOa_accept[iAJ] = new double[ndt];
        this->donOb_accept[iAJ] = new double[ndt];
        for (int it = 0; it < ndt; it++)
        {
            this->OOO_angle[iAJ][it] = 0;
            this->OH_central_angle[iAJ][it] = 0;
            this->OH_plane_angle[iAJ][it] = 0;
            this->OaOH_angle[iAJ][it] = 0;
            this->ObOH_angle[iAJ][it] = 0;
            this->OaO_distance[iAJ][it] = 0;
            this->ObO_distance[iAJ][it] = 0;
            this->central_anthH_donate[iAJ][it] = 0;
            this->central_accept[iAJ][it] = 0;
            this->donOa_donate[iAJ][it] = 0;
            this->donOb_donate[iAJ][it] = 0;
            this->donOa_accept[iAJ][it] = 0;
            this->donOb_accept[iAJ][it] = 0;
        }
    }
    //cout << "Allocation done." << endl;
    assert(INPUT.system == "water" or INPUT.system == "hydronium" or INPUT.system == "hydroxide");
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


void HBs_near_DoubleDonor::calc(CellFile &cel)
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
    //cout << "water setup done." << endl;
    for (int iDD=0; iDD<INPUT.nPT; iDD++)
    {
        if (this->DD[iDD].Oafter == -1 or this->DD[iDD].Obefore == -1)
        {
            continue;
        }
        if (abs(this->DD[iDD].DD_start_time - cel.snapshot_time) < INPUT.stay_tmax)
        {
            int ss_index = int((cel.snapshot_time - this->DD[iDD].DD_start_time)/INPUT.stay_dt);
            ss_index += int(this->ndt/2);
            assert(ss_index >= 0 and ss_index < this->ndt);
            if (iDD >= INPUT.nPT or ss_index >= this->ndt)
            {
                cout << iDD << " " << ss_index << endl;
                cout << "Warning! Array index out of range." << endl;
                exit(0);
            }
            int Ob = -1;
            if (this->DD[iDD].jump == "True")
            {
                Ob = DD[iDD].Oafter;
            }
            if (this->DD[iDD].jump == "False")
            {
                if (DD[iDD].Oab[0] == DD[iDD].Obefore)
                {
                    Ob = DD[iDD].Oab[1];
                }
                else
                {
                    Ob = DD[iDD].Oab[0];
                }
            }
            assert(Ob != -1 and Ob != DD[iDD].Obefore);
            this->central_accept[iDD][ss_index] = water[DD[iDD].centralO].naccept;
            int anthH_don = 0;
            for (int idon=0; idon<water[DD[iDD].centralO].ndonate; idon++)
            {
                if (water[DD[iDD].centralO].donateH[idon] != this->DD[iDD].H) anthH_don++;
            }
            this->central_anthH_donate[iDD][ss_index] = anthH_don;
            this->donOa_accept[iDD][ss_index] = water[DD[iDD].Obefore].naccept;
            this->donOa_donate[iDD][ss_index] = water[DD[iDD].Obefore].ndonate;
            this->donOb_accept[iDD][ss_index] = water[Ob].naccept;
            this->donOb_donate[iDD][ss_index] = water[Ob].ndonate;
            this->OaO_distance[iDD][ss_index] = distance(cel.atom[ito].pos[DD[iDD].Obefore], cel.atom[ito].pos[DD[iDD].centralO], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
            this->ObO_distance[iDD][ss_index] = distance(cel.atom[ito].pos[Ob], cel.atom[ito].pos[DD[iDD].centralO], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
            //cout << this->OaO_distance[iAJ][ss_index] << " " << this->ObO_distance[iAJ][ss_index] << endl;
            this->OaOH_angle[iDD][ss_index] = HBs::angle(cel, cel.atom[ito].pos[DD[iDD].Obefore], cel.atom[ito].pos[DD[iDD].centralO], cel.atom[ith].pos[DD[iDD].H]);
            //cout << "OaOH_angle[iAJ][ss_index] = " << OaOH_angle[iAJ][ss_index] << endl;
            this->ObOH_angle[iDD][ss_index] = HBs::angle(cel, cel.atom[ito].pos[Ob], cel.atom[ito].pos[DD[iDD].centralO], cel.atom[ith].pos[DD[iDD].H]);
            //cout << "ObOH_angle[iAJ][ss_index] = " << ObOH_angle[iAJ][ss_index] << endl;
            double central_angle = 0;
            double planar_angle = 0;
            HBs_near_AngularJump::calc_OH_central_plane_angle(cel, water, DD[iDD].centralO, DD[iDD].Obefore, Ob, DD[iDD].H, ito, ith, central_angle, planar_angle);
            //cout << "central_angle = " << central_angle << endl;
            this->OH_central_angle[iDD][ss_index] = central_angle;
            //cout << "OH_central_angle[iAJ][ss_index] = " << OH_central_angle[iAJ][ss_index] << endl;
            this->OH_plane_angle[iDD][ss_index] = planar_angle;
            //cout << "OH_plane_angle[iAJ][ss_index] = " << OH_plane_angle[iAJ][ss_index] << endl;
            //cout << AJ[iAJ].Oa << " " << centralO << " " << AJ[iAJ].Ob << endl;
            //cout << iAJ << " < " << INPUT.nPT << ", " << ss_index << " < " << this->ndt << endl;
            //cout << this->OOO_angle[iAJ][ss_index] << endl;
            this->OOO_angle[iDD][ss_index] = HBs::angle(cel, cel.atom[ito].pos[DD[iDD].Obefore], cel.atom[ito].pos[DD[iDD].centralO], cel.atom[ito].pos[Ob]);
            //cout << "OOO_angle[iAJ][ss_index] = " << OOO_angle[iAJ][ss_index] << endl;
            //cout << "Trial calculation done." << endl;
            //cout << "Over" << endl;
        }
    }
    delete[] water;
    return;
}

void HBs_near_DoubleDonor::output()
{
    this->single_output(central_accept, "central_accept_jump.txt", "central_accept_rattle.txt");
    this->single_output(central_anthH_donate, "central_anthH_donate_jump.txt", "central_anthH_donate_rattle.txt");
    this->single_output(donOa_donate, "donOa_donate_jump.txt", "donOa_donate_rattle.txt");
    this->single_output(donOb_donate, "donOb_donate_jump.txt", "donOb_donate_rattle.txt");
    this->single_output(donOa_accept, "donOa_accept_jump.txt", "donOa_accept_rattle.txt");
    this->single_output(donOb_accept, "donOb_accept_jump.txt", "donOb_accept_rattle.txt");
    this->single_output(OaO_distance, "OaO_distance_jump.txt", "OaO_distance_rattle.txt");
    this->single_output(ObO_distance, "ObO_distance_jump.txt", "ObO_distance_rattle.txt");
    this->single_output(OaOH_angle, "OaOH_angle_jump.txt", "OaOH_angle_rattle.txt");
    this->single_output(ObOH_angle, "ObOH_angle_jump.txt", "ObOH_angle_rattle.txt");
    this->single_output(OOO_angle, "OOO_angle_jump.txt", "OOO_angle_rattle.txt");
    this->single_output(OH_central_angle, "OH_central_angle_jump.txt", "OH_central_angle_rattle.txt");
    this->single_output(OH_plane_angle, "OH_plane_angle_jump.txt", "OH_plane_angle_rattle.txt");
    return;
}

void HBs_near_DoubleDonor::single_output(double** arr, string jump_file, string rattle_file)
{
    ofstream ofs_jump(jump_file);
    ofstream ofs_rattle(rattle_file);
    for (int it=0; it<this->ndt; it++)
    {
        double sum_jump = 0;
        double sum_rattle = 0;
        int njump = 0;
        int nrattle = 0;
        for (int iDD=0; iDD<INPUT.nPT; iDD++)
        {
            if (arr[iDD][it] != 0 and this->DD[iDD].jump == "True")
            {
                sum_jump += arr[iDD][it];
                njump++;
            }
            else if (arr[iDD][it] != 0 and this->DD[iDD].jump == "False")
            {
                sum_rattle += arr[iDD][it];
                nrattle++;
            }
        }
        if (njump > 0)
        {
            ofs_jump << it*INPUT.stay_dt - INPUT.stay_tmax << " " << sum_jump/njump << " " << njump << endl;
        }
        if (nrattle > 0)
        {
            ofs_rattle << it*INPUT.stay_dt - INPUT.stay_tmax << " " << sum_rattle/nrattle << " " << nrattle << endl;
        }
    }
    ofs_rattle.close();
    ofs_jump.close();
    return;
}