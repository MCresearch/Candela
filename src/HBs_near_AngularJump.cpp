#include "HBs_near_AngularJump.h"
#include "gfun.h"
HBs_near_AngularJump::HBs_near_AngularJump(){}

HBs_near_AngularJump::~HBs_near_AngularJump(){}

void HBs_near_AngularJump::Routine()
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
    AJ = new AngularJump[INPUT.nPT];
    AngularJump::setup_AJ(INPUT.nPT, this->AJ);
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
/*
void HBs_near_AngularJump::double_array_new(double** arr, int dim1, int dim2)
{
    arr = new double*[dim1];
    for (int i=0; i<dim1; i++)
    {
        arr[i] = new double[dim2];
        for (int j=0; j<dim2; j++)
        {
            arr[i][j] = 0;
        }
    }
    return;
}
*/
void HBs_near_AngularJump::calc(CellFile &cel)
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
    for (int iAJ=0; iAJ<INPUT.nPT; iAJ++)
    {
        int centralO = this->AJ[iAJ].centralO;
        //if (this->AJ[iAJ].Oa < 0 or this->AJ[iAJ].Ob < 0)
        //{
            if (this->AJ[iAJ].AJ_ss_time - cel.snapshot_time <= INPUT.relax_time and this->AJ[iAJ].AJ_ss_time - cel.snapshot_time >= INPUT.stay_tmax)
            {
                int OH_HB = 0;
                int Oa_index = -1;
                bool done = false;
                //cout << water[centralO].indexH[0] << " " << water[centralO].indexH[1] << " " << AJ[iAJ].H << endl;
                for (int idon=0; idon<water[centralO].ndonate; idon++)
                {
                    if(water[centralO].donateH[idon] == AJ[iAJ].H)
                    {
                        OH_HB++;
                        Oa_index = water[centralO].donateO[idon];
                        //cout << "Oa_index = " << Oa_index << ", OH_HB = " << OH_HB << endl;
                    }
                }
                if (OH_HB == 1)
                {
                    if (Oa_index == AJ[iAJ].Oab[0])
                    {
                        AJ[iAJ].Oa = Oa_index;
                        AJ[iAJ].Ob = AJ[iAJ].Oab[1];
                        done = true;
                    }
                    if (Oa_index == AJ[iAJ].Oab[1])
                    {
                        AJ[iAJ].Oa = Oa_index;
                        AJ[iAJ].Ob = AJ[iAJ].Oab[0];
                        done = true;
                    }
                }
            /*    if (OH_HB >= 2 or OH_HB == 0 or not done)
                {
                    //assert(Oa_index == AJ[iAJ].Oab[0] or Oa_index == AJ[iAJ].Oab[1]);
                    double dis0 = distance(cel.atom[ito].pos[centralO], cel.atom[ito].pos[AJ[iAJ].Oab[0]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                    double dis1 = distance(cel.atom[ito].pos[centralO], cel.atom[ito].pos[AJ[iAJ].Oab[1]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                    if (dis0 < dis1)
                    {
                        AJ[iAJ].Oa = AJ[iAJ].Oab[0];
                        AJ[iAJ].Ob = AJ[iAJ].Oab[1];
                    }
                    else
                    {
                        AJ[iAJ].Oa = AJ[iAJ].Oab[1];
                        AJ[iAJ].Ob = AJ[iAJ].Oab[0];
                    }
                }*/
            }
        //}

        if (abs(this->AJ[iAJ].AJ_ss_time - cel.snapshot_time) < INPUT.stay_tmax)
        {
            if (not ((AJ[iAJ].Oa == AJ[iAJ].Oab[0] and AJ[iAJ].Ob == AJ[iAJ].Oab[1]) or
            (AJ[iAJ].Oa == AJ[iAJ].Oab[1] and AJ[iAJ].Ob == AJ[iAJ].Oab[0])))
            {
                continue;
            }
            int ss_index = int((this->AJ[iAJ].AJ_ss_time - cel.snapshot_time)/INPUT.stay_dt);
            ss_index += int(this->ndt/2);
            assert(ss_index >= 0 and ss_index < this->ndt);
            if (iAJ >= INPUT.nPT or ss_index >= this->ndt)
            {
                cout << iAJ << " " << ss_index << endl;
                cout << "Warning! Array index out of range." << endl;
                exit(0);
            }
            this->central_accept[iAJ][ss_index] = water[centralO].naccept;
            int anthH_don = 0;
            for (int idon=0; idon<water[centralO].ndonate; idon++)
            {
                if (water[centralO].donateH[idon] != this->AJ[iAJ].H) anthH_don++;
            }
            this->central_anthH_donate[iAJ][ss_index] = anthH_don;
            this->donOa_accept[iAJ][ss_index] = water[AJ[iAJ].Oa].naccept;
            this->donOa_donate[iAJ][ss_index] = water[AJ[iAJ].Oa].ndonate;
            this->donOb_accept[iAJ][ss_index] = water[AJ[iAJ].Ob].naccept;
            this->donOb_donate[iAJ][ss_index] = water[AJ[iAJ].Ob].ndonate;
            //cout << "this->donOb_donate[iAJ][ss_index] = " << this->donOb_donate[iAJ][ss_index] << endl;
            this->OaO_distance[iAJ][ss_index] = distance(cel.atom[ito].pos[AJ[iAJ].Oa], cel.atom[ito].pos[centralO], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
            this->ObO_distance[iAJ][ss_index] = distance(cel.atom[ito].pos[AJ[iAJ].Ob], cel.atom[ito].pos[centralO], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
            //cout << this->OaO_distance[iAJ][ss_index] << " " << this->ObO_distance[iAJ][ss_index] << endl;
            this->OaOH_angle[iAJ][ss_index] = HBs::angle(cel, cel.atom[ito].pos[AJ[iAJ].Oa], cel.atom[ito].pos[centralO], cel.atom[ith].pos[AJ[iAJ].H]);
            //cout << "OaOH_angle[iAJ][ss_index] = " << OaOH_angle[iAJ][ss_index] << endl;
            this->ObOH_angle[iAJ][ss_index] = HBs::angle(cel, cel.atom[ito].pos[AJ[iAJ].Ob], cel.atom[ito].pos[centralO], cel.atom[ith].pos[AJ[iAJ].H]);
            //cout << "ObOH_angle[iAJ][ss_index] = " << ObOH_angle[iAJ][ss_index] << endl;
            double central_angle = 0;
            double planar_angle = 0;
            calc_OH_central_plane_angle(cel, water, centralO, AJ[iAJ].Oa, AJ[iAJ].Ob, AJ[iAJ].H, ito, ith, central_angle, planar_angle);
            //cout << "central_angle = " << central_angle << endl;
            this->OH_central_angle[iAJ][ss_index] = central_angle;
            //cout << "OH_central_angle[iAJ][ss_index] = " << OH_central_angle[iAJ][ss_index] << endl;
            this->OH_plane_angle[iAJ][ss_index] = planar_angle;
            //cout << "OH_plane_angle[iAJ][ss_index] = " << OH_plane_angle[iAJ][ss_index] << endl;
            //cout << AJ[iAJ].Oa << " " << centralO << " " << AJ[iAJ].Ob << endl;
            //cout << iAJ << " < " << INPUT.nPT << ", " << ss_index << " < " << this->ndt << endl;
            //cout << this->OOO_angle[iAJ][ss_index] << endl;
            this->OOO_angle[iAJ][ss_index] = HBs::angle(cel, cel.atom[ito].pos[AJ[iAJ].Oa], cel.atom[ito].pos[centralO], cel.atom[ito].pos[AJ[iAJ].Ob]);
            //cout << "OOO_angle[iAJ][ss_index] = " << OOO_angle[iAJ][ss_index] << endl;
            //cout << "Trial calculation done." << endl;
            //cout << "Over" << endl;
        }
    }
    delete[] water;
    return;
}


void HBs_near_AngularJump::calc_OH_central_plane_angle(Cell &cel, Water* water, const int &iwater, const int &iwater1, const int &iwater2,  const int &iH, const int &ito, const int &ith, double &central_angle, double &planar_angle)
{
    Vector3<double> posO1;
    Vector3<double> posO2;
    Vector3<double> posH;
    posO1 = cel.atom[ito].pos[iwater1];
    posO2 = cel.atom[ito].pos[iwater2];
    posH = cel.atom[ith].pos[iH];
    Dist2::putback_cell(cel.atom[ito].pos[iwater], posO1);
    Dist2::putback_cell(cel.atom[ito].pos[iwater], posO2);
    Dist2::putback_cell(cel.atom[ito].pos[iwater], posH);
    Vector3<double> vec_a = posO1 - cel.atom[ito].pos[iwater];
    Vector3<double> vec_b = posO2 - cel.atom[ito].pos[iwater];
    Vector3<double> vec_x = posH - cel.atom[ito].pos[iwater];
    vec_a = vec_a/vec_a.norm();
    vec_b = vec_b/vec_b.norm();
    vec_x = vec_x/vec_x.norm();
    double ab_prod = dot(vec_a, vec_b);
    double ax_prod = dot(vec_a, vec_x);
    double bx_prod = dot(vec_b, vec_x);

    double A = (ax_prod + bx_prod)/(1 + ab_prod)/2 + (ax_prod - bx_prod)/(1-ab_prod)/2;
    double B = (ax_prod + bx_prod)/(1 + ab_prod)/2 - (ax_prod - bx_prod)/(1-ab_prod)/2;

    Vector3<double> proj_x = A*vec_a + B*vec_b;
    Vector3<double> central_arrow = (vec_a+vec_b)/(vec_a+vec_b).norm();
    Vector3<double> vertical_x = vec_x - proj_x;
    //cout << dot(proj_x, vertical_x) << endl;
    //cout << "proj = " << proj_x.x << " " << proj_x.y << " " << proj_x.z << endl;
    //cout << "vertical = " << vertical_x.x << " " << vertical_x.y << " " << vertical_x.z << endl;
    assert(abs(dot(proj_x, vertical_x)) < 1e-5);
    planar_angle = dot(proj_x, vec_x)/vec_x.norm()/proj_x.norm();
    planar_angle = acos(planar_angle) / PI * 180;
    central_angle = dot(central_arrow, proj_x)/proj_x.norm();
    central_angle = acos(central_angle) / PI * 180;
    double proj_a_angle = dot(proj_x, vec_a)/proj_x.norm()/vec_a.norm();
    proj_a_angle = acos(proj_a_angle) / PI*180;
    double proj_b_angle = dot(proj_x, vec_b)/proj_x.norm()/vec_b.norm();
    proj_b_angle = acos(proj_b_angle) / PI*180;
    if (proj_a_angle > proj_b_angle) central_angle = -central_angle;
    return;
}

void HBs_near_AngularJump::output()
{
    single_output(this->central_accept, INPUT.nPT, this->ndt, "central_accept.txt");
    single_output(this->central_anthH_donate, INPUT.nPT, this->ndt, "central_anthH_donate.txt");
    single_output(this->donOa_accept, INPUT.nPT, this->ndt, "donOa_accept.txt");
    single_output(this->donOa_donate, INPUT.nPT, this->ndt, "donOa_donate.txt");
    single_output(this->donOb_accept, INPUT.nPT, this->ndt, "donOb_accept.txt");
    single_output(this->donOb_donate, INPUT.nPT, this->ndt, "donOb_donate.txt");
    single_output(this->OaO_distance, INPUT.nPT, this->ndt, "OaO_distance.txt");
    single_output(this->ObO_distance, INPUT.nPT, this->ndt, "ObO_distance.txt");
    single_output(this->OaOH_angle, INPUT.nPT, this->ndt, "OaOH_angle.txt");
    single_output(this->ObOH_angle, INPUT.nPT, this->ndt, "ObOH_angle.txt");
    single_output(this->OH_central_angle, INPUT.nPT, this->ndt, "OH_central_angle.txt");
    single_output(this->OH_plane_angle, INPUT.nPT, this->ndt, "OH_planar_angle.txt");
    single_output(this->OOO_angle, INPUT.nPT, this->ndt, "OOO_angle.txt");
    return;
}

void HBs_near_AngularJump::single_output(double** arr, int dim1, int dim2, string file_name)
{
    ofstream ofs(file_name);
    for (int it=0; it<dim2; it++)
    {
        double sum = 0;
        int nvalid = 0;
        for (int ii = 0; ii<dim1; ii++)
        {
            if (arr[ii][it] != 0)
            {
                nvalid++;
                sum += arr[ii][it];
            }
        }
        if (nvalid > 0)
        {
            sum /= nvalid;
            ofs << it*INPUT.stay_dt - INPUT.stay_tmax << " " << sum << endl;
        }
    }
    ofs.close();
    return;
}