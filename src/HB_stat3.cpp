#include "HB_stat3.h"
#include "dist2.h"
HB_stat3::HB_stat3(){}
HB_stat3::~HB_stat3(){}

void HB_stat3::Routine()
{
    TITLE("HB_stat3","Routine");
    ofs_running << "HB reorientation statistics." << endl;
    ofs_running << "func = 1: search for all the reorientation incidents." << endl;
    ofs_running << "func = 2: analyze the chrononical properties of reorientation." << endl;

    if (INPUT.func == 1)
    {
        int nangle = int(INPUT.theta/INPUT.dtheta);
        this->OOH_angle = new double[nangle];
        this->OOO_angle = new double[nangle];
        this->OH_central_angle = new double[nangle];
        this->OH_plane_angle = new double[nangle];
        for (int ia = 0; ia<nangle; ia++)
        {
            OOH_angle[ia] = 0;
            OOO_angle[ia] = 0;
            OH_central_angle[ia] = 0;
            OH_plane_angle[ia] = 0;
        }
        this->OO_distance = new double[int(INPUT.rcut/INPUT.dr)];
        for (int ia = 0; ia<int(INPUT.rcut/INPUT.dr); ia++)
        {
            OO_distance[ia] = 0;
        }
        this->central_accept = new double[10];
        this->central_donate = new double[10];  
        this->donO_accept = new double[10];
        this->donO_donate = new double[10];
        for (int idon=0; idon<10; idon++)
        {
            central_accept[idon] = 0;
            central_donate[idon] = 0;
            donO_accept[idon] = 0;
            donO_donate[idon] = 0;
        }
        this->ndouble_donate = 0;
        this->don_acc = new double*[10];
        for (int idon=0; idon<10; idon++)
        {
            this->don_acc[idon] = new double[10];
            for (int iacc=0; iacc<10; iacc++)
            {
                this->don_acc[idon][iacc] = 0;
            }
        }
        this->ofs_mid.open("middle_log.txt");
    }


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
    ofs_mid.close();
    this->normalize_output(count_geometry_number);
    return;
}

void HB_stat3::calc(Cell &cel)
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
        if (water[iwater].ndonate < 10 and water[iwater].naccept < 10)
        {
            don_acc[water[iwater].ndonate][water[iwater].naccept]++;
        }
        if (water[iwater].ndonate == 1)
        {
            continue;
        }
        for (int idon1=0; idon1<water[iwater].ndonate; idon1++)
        {
            for (int idon2=idon1+1; idon2<water[iwater].ndonate; idon2++)
            {
                if (water[iwater].donateH[idon1] != water[iwater].donateH[idon2]) continue;
                int iwater1 = water[iwater].donateO[idon1];
                int iwater2 = water[iwater].donateO[idon2];
                int iH = water[iwater].donateH[idon2];
                this->ndouble_donate++;
                double OOO = HBs::angle(cel, cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2]);
                if (OOO < INPUT.theta)
                {
                    this->OOO_angle[int(OOO/INPUT.dtheta)]++;
                }
                double OOH1 = HBs::angle(cel, cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater], cel.atom[ith].pos[iH]);
                if (OOH1 < INPUT.theta)
                {
                    this->OOH_angle[int(OOH1/INPUT.dtheta)]++;
                }
                double OOH2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater], cel.atom[ith].pos[iH]);
                if (OOH2 < INPUT.theta)
                {
                    this->OOH_angle[int(OOH2/INPUT.dtheta)]++;
                }
                double central_angle = this->calc_OH_central_plane_angle(cel, water, iwater, iwater1, iwater2, iH, ito, ith);
                double OO_distance1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater1], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                double OO_distance2 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                if (OO_distance1 < INPUT.rcut)
                {
                    this->OO_distance[int(OO_distance1/INPUT.dr)]++;
                }
                if (OO_distance2 < INPUT.rcut)
                {
                    this->OO_distance[int(OO_distance2/INPUT.dr)]++;
                }

                this->central_accept[water[iwater].naccept]++;
                this->central_donate[water[iwater].ndonate]++;
                this->donO_accept[water[iwater1].naccept]++;
                this->donO_accept[water[iwater2].naccept]++;
                this->donO_donate[water[iwater1].ndonate]++;
                this->donO_donate[water[iwater2].ndonate]++;
                if (central_angle < 3)
                {
                    this->ofs_mid << setprecision(10) << setw(20) << cel.snapshot_index << setw(20) << cel.snapshot_time << setw(20)
                     << setw(20) << iwater;
                    if (distance(cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3) < 
                    distance(cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3))
                    {
                        this->ofs_mid << setw(20) << iwater1 << setw(20) << iwater2;
                    }
                    else
                    {
                        this->ofs_mid << setw(20) << iwater2 << setw(20) << iwater1;
                    }
                    this->ofs_mid << setw(20) << central_angle << setw(20) << iH <<endl;    
                }
            }
        }
    }

    delete[] water;
}

double HB_stat3::calc_OH_central_plane_angle(Cell &cel, Water* water, const int &iwater, const int &iwater1, const int &iwater2,  const int &iH, const int &ito, const int &ith)
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
    double plane_angle = dot(proj_x, vec_x)/vec_x.norm()/proj_x.norm();
    plane_angle = acos(plane_angle) / PI * 180;
    if (plane_angle < INPUT.theta)
    {
        this->OH_plane_angle[int(plane_angle/INPUT.dtheta)]++;
    }
    double central_angle = dot(central_arrow, proj_x)/proj_x.norm();
    central_angle = acos(central_angle) / PI * 180;
    if (central_angle < INPUT.theta)
    {
        this->OH_central_angle[int(central_angle/INPUT.dtheta)]++;
    }
    return central_angle;
}

void HB_stat3::normalize_output(int &count_geometry_number)
{
    double sum = 0;
    for (int idon=0; idon<10; idon++)
    {
        for(int iacc=0; iacc<10; iacc++)
        {
            sum += this->don_acc[idon][iacc] ;
        }
    }
    for (int idon=0; idon<10; idon++)
    {
        for(int iacc=0; iacc<10; iacc++)
        {
            this->don_acc[idon][iacc] /= sum;
        }
    }
    for (int idon=0; idon<10; idon++)
    {
        for(int iacc=0; iacc<10; iacc++)
        {
            ofs_running << setw(20) << don_acc[idon][iacc];
        }
        ofs_running << endl;
    }
    single_normalize(this->OOO_angle, INPUT.theta, INPUT.dtheta);
    single_normalize(this->OOH_angle, INPUT.theta, INPUT.dtheta);
    single_normalize(this->OH_central_angle, INPUT.theta, INPUT.dtheta);
    single_normalize(this->OH_plane_angle, INPUT.theta, INPUT.dtheta);
    single_normalize(this->OO_distance, INPUT.rcut, INPUT.dr);
    double maxn = 10;
    double one = 1;
    single_normalize(this->central_accept, maxn, one);
    single_normalize(this->central_donate, maxn, one);
    single_normalize(this->donO_donate, maxn, one);
    single_normalize(this->donO_accept, maxn, one);

    this->single_output(this->OOO_angle, INPUT.theta, INPUT.dtheta, "OOO_angle.txt", true);
    this->single_output(this->OOH_angle, INPUT.theta, INPUT.dtheta, "OOH_angle.txt", true);
    this->single_output(this->OH_central_angle, INPUT.theta, INPUT.dtheta, "central_angle.txt", true);
    this->single_output(this->OH_plane_angle, INPUT.theta, INPUT.dtheta, "plane_angle.txt", true);
    this->single_output(this->OO_distance, INPUT.rcut, INPUT.dr, "OO_distance.txt", true);
    this->single_output(this->central_accept, maxn, one, "central_accept.txt", false);
    this->single_output(this->central_donate, maxn, one, "central_donate.txt", false);
    this->single_output(this->donO_accept, maxn, one, "donO_accept.txt", false);
    this->single_output(this->donO_donate, maxn, one, "donO_donate.txt", false);
    ofs_running << "Number of double H donor = " << this->ndouble_donate << endl;
    return;
}

void HB_stat3::single_normalize(double* arr, double &rcut, double &dr)
{
    double sum_=0;
    for (int ir=0; ir<int(rcut/dr); ir++)
    {
        sum_ += arr[ir];
    }
    sum_ *= dr;
    for (int ir=0; ir<int(rcut/dr); ir++)
    {
        arr[ir] /= sum_;
    }
    return;
}

void HB_stat3::single_output(double* arr, double &rcut, double &dr, string file_name, bool half_plus)
{
    ofstream ofs_file(file_name);
    for (int ir=0; ir < int(rcut/dr); ir++)
    {
        if (half_plus)
        {
            ofs_file << (ir+0.5)*dr << " " << arr[ir] << endl;
        }
        else
        {
            ofs_file << ir*dr << " " << arr[ir] << endl;
        }
    }
    ofs_file.close();
    return;
}