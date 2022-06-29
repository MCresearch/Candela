#include "first_shell_angle.h"

First_shell_angle::First_shell_angle()
{}

First_shell_angle::~First_shell_angle()
{}

void First_shell_angle::Routine()
{
    this->nangle = int(180/INPUT.bdf_dtheta)+1;
    this->A_OOc_angle = new double[nangle];
    this->D_HOO_angle = new double[nangle];
    this->A_alpha_angle = new double[nangle];
    this->A_beta_angle = new double[nangle];
    this->D_alpha_angle = new double[nangle];
    this->D_beta_angle = new double[nangle];

    for (int iangle=0; iangle<this->nangle; iangle++)
    {
        this->A_OOc_angle[iangle] = 0;
        this->A_alpha_angle[iangle] = 0;
        this->A_beta_angle[iangle] = 0;
        this->D_HOO_angle[iangle] = 0;
        this->D_alpha_angle[iangle] = 0;
        this->D_beta_angle[iangle] = 0;
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
    this->output();
    return;
}

void First_shell_angle::calc(CellFile &cel)
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
        if (water[iwater].ndonate != INPUT.ndon or water[iwater].naccept != INPUT.nacc) continue;

        Vector3<double> pos_H1 = cel.atom[ith].pos[water[iwater].indexH[0]];
        Vector3<double> pos_H2 = cel.atom[ith].pos[water[iwater].indexH[1]];
        Dist2::putback_cell(cel.atom[ito].pos[iwater], pos_H1);
        Dist2::putback_cell(cel.atom[ito].pos[iwater], pos_H2);

        /*
        cout << "iwater = " << iwater << endl;
        cout << "O1 pos = " << cel.atom[ito].pos[iwater].x << " " << cel.atom[ito].pos[iwater].y << " " << cel.atom[ito].pos[iwater].z << endl;
        cout << "H1 pos = " << pos_H1.x << " " << pos_H1.y << " " << pos_H1.z << endl;
        cout << "H2 pos = " << pos_H2.x << " " << pos_H2.y << " " << pos_H2.z << endl;
        */
        pos_H1 = pos_H1 - cel.atom[ito].pos[iwater];
        pos_H2 = pos_H2 - cel.atom[ito].pos[iwater];
        pos_H1 = pos_H1/pos_H1.norm();
        pos_H2 = pos_H2/pos_H2.norm();
        Vector3<double> pos_center = pos_H1 + pos_H2;

        Vector3<double> vec_center_o = pos_center/pos_center.norm();
        //cout << "vec_center_o = " << vec_center_o.x << " " << vec_center_o.y << " " << vec_center_o.z << endl;
        Vector3<double> vec_per = cross(pos_H1, pos_H2);
        vec_per = vec_per/vec_per.norm();
        // vec_center_o and vec_per are centered on O atom.
        Vector3<double> vec_par = cross(vec_per, vec_center_o);
        vec_par = vec_par/vec_par.norm();

        Vector3<double> pos_center_up = cel.atom[ito].pos[iwater] + pos_center;
        pos_center = cel.atom[ito].pos[iwater] - pos_center;

        for (int idon=0; idon<water[iwater].ndonate; idon++)
        {
            int iwater2 = water[iwater].donateO[idon];
            int hindex_don = water[iwater].donateH[idon];
            double angle_hoo = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater], cel.atom[ith].pos[hindex_don]);
            this->D_HOO_angle[int(angle_hoo/INPUT.bdf_dtheta)]++;

            Vector3<double> vec_O2 = cel.atom[ito].pos[iwater2];
            Dist2::putback_cell(cel.atom[ito].pos[iwater], vec_O2);
            vec_O2 = vec_O2 - cel.atom[ito].pos[iwater];
            double x_per = vec_O2 * vec_per;
            double y_cen = vec_O2 * vec_center_o;
            Vector3<double> pos_O2_proj = cel.atom[ito].pos[iwater] + x_per*vec_per + y_cen*vec_center_o;
            double angle_alpha = HBs::angle(cel, pos_center_up, cel.atom[ito].pos[iwater], pos_O2_proj);
            this->D_alpha_angle[int(angle_alpha/INPUT.bdf_dtheta)]++;

            double x_par = vec_O2 * vec_par;
            Vector3<double> pos_O2_proj2 = cel.atom[ito].pos[iwater] + x_par*vec_par + y_cen*vec_center_o;
            double angle_beta = HBs::angle(cel, pos_center_up, cel.atom[ito].pos[iwater], pos_O2_proj2);
            this->D_beta_angle[int(angle_beta/INPUT.bdf_dtheta)]++;
        }

        for (int iacc=0; iacc<water[iwater].naccept; iacc++)
        {
            int iwater2 = water[iwater].acceptO[iacc];
            //cout << "pos of iwater2 " << iwater2 << " = " << cel.atom[ito].pos[iwater2].x << " " << cel.atom[ito].pos[iwater2].y << " " << cel.atom[ito].pos[iwater2].z << endl;
            double angle_center = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater], pos_center);
            this->A_OOc_angle[int(angle_center/INPUT.bdf_dtheta)]++;

            Vector3<double> vec_O2 = cel.atom[ito].pos[iwater2];
            Dist2::putback_cell(cel.atom[ito].pos[iwater], vec_O2);
            vec_O2 = vec_O2 - cel.atom[ito].pos[iwater];
            double x_per = vec_O2 * vec_per;
            double y_cen = vec_O2 * vec_center_o;
            Vector3<double> pos_O2_proj = cel.atom[ito].pos[iwater] + x_per*vec_per + y_cen*vec_center_o;
            double angle_alpha = HBs::angle(cel, pos_center, cel.atom[ito].pos[iwater], pos_O2_proj);
            this->A_alpha_angle[int(angle_alpha/INPUT.bdf_dtheta)]++;

            double x_par = vec_O2 * vec_par;
            Vector3<double> pos_O2_proj2 = cel.atom[ito].pos[iwater] + x_par*vec_par + y_cen*vec_center_o;
            double angle_beta = HBs::angle(cel, pos_center, cel.atom[ito].pos[iwater], pos_O2_proj2);
            this->A_beta_angle[int(angle_beta/INPUT.bdf_dtheta)]++;
        }
    }

    delete[] water;
}

void First_shell_angle::output()
{
    double summ_D = 0;
    double summ_A = 0;
    double summ_A_a = 0;
    double summ_A_b = 0;
    double summ_D_a = 0;
    double summ_D_b = 0;
    for (int iangle=0; iangle<int(180/INPUT.bdf_dtheta)+1; iangle++)
    {
        summ_D += this->D_HOO_angle[iangle];
        summ_A += this->A_OOc_angle[iangle];
        summ_A_a += this->A_alpha_angle[iangle];
        summ_A_b += this->A_beta_angle[iangle];
        summ_D_a += this->D_alpha_angle[iangle];
        summ_D_b += this->D_beta_angle[iangle];
    }

    ofstream ofs("D_HOO_angle.txt");
    ofstream ofs2("A_OOc_angle.txt");
    ofstream ofs3("A_alpha_angle.txt");
    ofstream ofs4("A_beta_angle.txt");
    ofstream ofs5("D_alpha_angle.txt");
    ofstream ofs6("D_beta_angle.txt");
    for (int iangle=0; iangle<int(180/INPUT.bdf_dtheta)+1; iangle++)
    {
        ofs << iangle*INPUT.bdf_dtheta << " " << this->D_HOO_angle[iangle]/summ_D/INPUT.bdf_dtheta << endl;
        ofs2 << iangle*INPUT.bdf_dtheta << " " << this->A_OOc_angle[iangle]/summ_A/INPUT.bdf_dtheta << endl;
        ofs3 << iangle*INPUT.bdf_dtheta << " " << this->A_alpha_angle[iangle]/summ_A_a/INPUT.bdf_dtheta << endl;
        ofs4 << iangle*INPUT.bdf_dtheta << " " << this->A_beta_angle[iangle]/summ_A_b/INPUT.bdf_dtheta << endl;
        ofs5 << iangle*INPUT.bdf_dtheta << " " << this->D_alpha_angle[iangle]/summ_D_a/INPUT.bdf_dtheta << endl;
        ofs6 << iangle*INPUT.bdf_dtheta << " " << this->D_beta_angle[iangle]/summ_D_b/INPUT.bdf_dtheta << endl;

    }
    ofs.close();
    ofs2.close();
    ofs3.close();
    ofs4.close();
    ofs5.close();
    ofs6.close();
    return;
}