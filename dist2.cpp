#include "dist2.h"
#include "cellFile.h"
#include "input.h"
#include "HBs.h"
#include "Wan_centers_stat.h"
Dist2::Dist2(){}

Dist2::~Dist2(){}

void Dist2::Routine()
{
	cout << "Compute 3D distribution of atoms/MLWFs around water molecule." << endl;
	cout << "The format will be charge density format or wave function format that can be read by VESTA." << endl;
    cout << "The main purpose of this code is to show the MLWF shape change and 2nd solvation shell shape change of water under pressure." << endl;
    
    assert(INPUT.u1>0);

    this->dr = INPUT.rcut / INPUT.u1;

    cout << "dr = " << dr << endl;

    this->dist3D = new double** [INPUT.u1];
    this->dist3D2 = new double** [INPUT.u1];
    for (int iu=0; iu<INPUT.u1; iu++)
    {
        this->dist3D[iu] = new double*[INPUT.u1];
        this->dist3D2[iu] = new double*[INPUT.u1];
        for (int iu1=0; iu1<INPUT.u1; iu1++)
        {
            this->dist3D[iu][iu1] = new double[INPUT.u1];
            this->dist3D2[iu][iu1] = new double[INPUT.u1];
            for (int iu2=0; iu2<INPUT.u1; iu2++)
            {
                this->dist3D[iu][iu1][iu2] = 0.0;
                this->dist3D2[iu][iu1][iu2] = 0.0;
            }//iu2
        }//iu1
    }//iu

    this->which_HOO_angle = new double[6]; 
    // 0: thisHD->D; 1: otherHD->D; 2: A->D; 3: thisHD->A; 4: otherHD->A; 5: A->A.
    for (int iwhich=0; iwhich<6; iwhich++)
    {
        which_HOO_angle[iwhich] = 0;
    }
	this->count_geometry_number=0;
    this->nabnormal = 0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo)
	{
		//cout << "count_geometry_number = " << count_geometry_number << endl;
        // cel : input geometry file
		CellFile cel;
        cout << "igeo=" << igeo << endl;
		if(igeo%INPUT.geo_interval!=0 or igeo < INPUT.geo_ignore) cel.read_and_used=false;
		else cel.read_and_used=true;

		stringstream ss; ss << igeo;
		cel.file_name = ss.str();

		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) 
        {
            cel.clean();
            continue;
        }

        this->calculate_3D(cel);
        cout << "count_geometry_number = " << count_geometry_number << endl;
        cel.clean();
        count_geometry_number++;
    }

    double sum = 0;
    for(int ix=0; ix<INPUT.u1; ix++)
    {
        for(int iy=0; iy<INPUT.u1; iy++)
        {
            for(int iz=0; iz<INPUT.u1; iz++)
            {
                sum += dist3D[ix][iy][iz];
            }
        }
    }
    cout << "sum = " << sum << endl;
    assert(sum > 0);
    for(int iz=0; iz<INPUT.u1; iz++)
    {
        for(int iy=0; iy<INPUT.u1; iy++)
        {
            for(int ix=0; ix<INPUT.u1; ix++)
            {
                dist3D[ix][iy][iz] /= sum;
                //cout << dist3D[ix][iy][iz] << endl;
            }
        }
    }

    double sum2=0;
    if (INPUT.func == 2)
    {
        for(int ix=0; ix<INPUT.u1; ix++)
        {
            for(int iy=0; iy<INPUT.u1; iy++)
            {
                for(int iz=0; iz<INPUT.u1; iz++)
                {
                    sum2 += dist3D2[ix][iy][iz];
                }
            }
        }
        cout << "sum2 = " << sum2 << endl;
        assert(sum2 > 0);
        for(int ix=0; ix<INPUT.u1; ix++)
        {
            for(int iy=0; iy<INPUT.u1; iy++)
            {
                for(int iz=0; iz<INPUT.u1; iz++)
                {
                    dist3D2[ix][iy][iz] /= sum2;
                }
            }
        }
    }

    ofstream ofs("dist3D.xsf");
    ofstream ofs_txt("dist3D.txt");
    ofs_txt << "dist_3D" << endl;
    ofs << "CRYSTAL\nPRIMVEC\n";
    ofs << INPUT.rcut << " " << 0 << " " << 0 << endl;
    ofs << 0 << " " << INPUT.rcut << " " << 0 << endl;
    ofs << "0 0 " << INPUT.rcut << endl;
    ofs << "PRIMCOORD" << endl;
    ofs << "3 1" << endl;
    ofs << "O " << INPUT.rcut/2 << " " << INPUT.rcut/2 << " " << INPUT.rcut/2 << endl;
    if (INPUT.ndon != 1)
    {
        ofs << "H " << INPUT.rcut/2+sin(PI*52.5/180) << " " << INPUT.rcut/2+cos(PI*52.5/180) << " " << INPUT.rcut/2 << endl;
        ofs << "H " << INPUT.rcut/2-sin(PI*52.5/180) << " " << INPUT.rcut/2+cos(PI*52.5/180) << " " << INPUT.rcut/2 << endl;
    }
    else
    {
        ofs << "H " << INPUT.rcut/2+1 << " " << INPUT.rcut/2 << " " << INPUT.rcut/2 << endl;
        ofs << "H " << INPUT.rcut/2-sin(PI*15/180) << " " << INPUT.rcut/2+cos(PI*15/180) << " " << INPUT.rcut/2 << endl;
    }
    ofs << " BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n DATAGRID_3D_UNKNOWN\n";
    ofs << " " << INPUT.u1 << " " << INPUT.u1 << " " << INPUT.u1 << endl;
    ofs << " 0 0 0" << endl;
    ofs << INPUT.rcut << " " << 0 << " " << 0 << endl;
    ofs << 0 << " " << INPUT.rcut << " " << 0 << endl;
    ofs << "0 0 " << INPUT.rcut << endl;
    for(int iz=0; iz<INPUT.u1; ++iz)
	{
		for(int iy=0; iy<INPUT.u1; ++iy)
		{
			for(int ix=0; ix<INPUT.u1; ++ix)
			{
				ofs << dist3D[ix][iy][iz] << endl;
                ofs_txt << dist3D[ix][iy][iz] << endl;
			}
		}
	}
	ofs << " END_DATAGRID_3D" << endl;
	ofs << " END_BLOCK_DATAGRID_3D" << endl;
	ofs.close();
    ofs_txt.close();
    if (INPUT.func == 2)
    {
        ofstream ofs2("dist3D2.xsf");
        ofstream ofs2_txt("dist3D2.txt");
        ofs2_txt << "dist_3D" << endl;
        ofs2 << "CRYSTAL\nPRIMVEC\n";
        ofs2 << INPUT.rcut << " " << 0 << " " << 0 << endl;
        ofs2 << 0 << " " << INPUT.rcut << " " << 0 << endl;
        ofs2 << "0 0 " << INPUT.rcut << endl;
        ofs2 << "PRIMCOORD" << endl;
        ofs2 << "3 1" << endl;
        ofs2 << "O " << INPUT.rcut/2 << " " << INPUT.rcut/2 << " " << INPUT.rcut/2 << endl;
        ofs2 << "H " << INPUT.rcut/2+1 << " " << INPUT.rcut/2 << " " << INPUT.rcut/2 << endl;
        ofs2 << "H " << INPUT.rcut/2-sin(PI*15/180) << " " << INPUT.rcut/2+cos(PI*15/180) << " " << INPUT.rcut/2 << endl;
        ofs2 << " BEGIN_BLOCK_DATAGRID_3D\n 3D_PWSCF\n DATAGRID_3D_UNKNOWN\n";
        ofs2 << " " << INPUT.u1 << " " << INPUT.u1 << " " << INPUT.u1 << endl;
        ofs2 << " 0 0 0" << endl;
        ofs2 << INPUT.rcut << " " << 0 << " " << 0 << endl;
        ofs2 << 0 << " " << INPUT.rcut << " " << 0 << endl;
        ofs2 << "0 0 " << INPUT.rcut << endl;
        for(int iz=0; iz<INPUT.u1; ++iz)
        {
            for(int iy=0; iy<INPUT.u1; ++iy)
            {
                for(int ix=0; ix<INPUT.u1; ++ix)
                {
                    ofs2 << dist3D2[ix][iy][iz] << endl;
                    ofs2_txt << dist3D2[ix][iy][iz] << endl;
                }
            }
        }
        ofs2 << " END_DATAGRID_3D" << endl;
        ofs2 << " END_BLOCK_DATAGRID_3D" << endl;
        ofs2.close();
        ofs2_txt.close();
    }

    if (INPUT.func == 2)
    {
        sum = 0;
        for (int i=0; i<6; i++)
        {
            sum += this->which_HOO_angle[i];
        }
        assert(sum > 0);
        for (int i=0; i<6; i++)
        {
            this->which_HOO_angle[i] /= sum;
        }
        ofs_running << "The distribution of H-O-O angle types:" << endl;
        ofs_running << setw(20) << "first_shell " << setw(20) << "second_shell " << setw(20) << "percentage (%)" << endl;
        ofs_running << setw(20) << "this_H_donate" << setw(20) << "donate" << setw(20) << 100*this->which_HOO_angle[0] << endl;
        ofs_running << setw(20) << "another_H_donate" << setw(20) << "donate" << setw(20) << 100*this->which_HOO_angle[1] << endl;
        ofs_running << setw(20) << "accept" << setw(20) << "donate" << setw(20) << 100*this->which_HOO_angle[2] << endl;
        ofs_running << setw(20) << "this_H_donate" << setw(20) << "accept" << setw(20) << 100*this->which_HOO_angle[3] << endl;
        ofs_running << setw(20) << "another_H_donate" << setw(20) << "accept" << setw(20) << 100*this->which_HOO_angle[4] << endl;
        ofs_running << setw(20) << "accept" << setw(20) << "accept" << setw(20) << 100*this->which_HOO_angle[5] << endl;
    }
    

    cout << "number of water molecules with acute angle = " << this->nabnormal << endl;

    return;

}

void Dist2::calculate_3D(Cell &cel)
{
    assert(INPUT.system == "water");

	int ito=-1;
	int ith=-1;
	int itcl=-1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
		else if(cel.atom[it].id=="Cl") itcl=it;
	}
	if(INPUT.ntype>=2){ assert(ito>=0); assert(ith>=0);}

	Water *water = new Water[cel.atom[ito].na];
	Water::nions=0;

    int** bond_wan_index = new int*[cel.atom[ito].na];
    int** lone_wan_index = new int*[cel.atom[ito].na];
    for (int ia=0; ia<cel.atom[ito].na; ia++)
    {
        bond_wan_index[ia] = new int[3];
        lone_wan_index[ia] = new int[3];
    }

	HBs::setup_water(cel, water);
    if (INPUT.func == 1)// func = 1: calculate for Wannier functions.
    {
        Wan_centers_stat::allocate_wan(water, cel, bond_wan_index, lone_wan_index, ito, ith);
    }
    for (int iwater=0; iwater<cel.atom[ito].na; iwater++)
    {
        // for func = 2
        int* first_shell_index = new int[6];
        int** second_shell_index = new int*[6];
        int* nwater_second_shell = new int[6];
        for (int ish=0; ish<6; ish++)
        {
            first_shell_index[ish] = -1;
            second_shell_index[ish] = new int[6];
            nwater_second_shell[ish] = 0;
            for (int ish2=0; ish2<6; ish2++)
            {
                second_shell_index[ish][ish2] = -1;
            }
        }
        int nwater_first_shell = 0;

        if (INPUT.func == 1)
        {
            Vector3<double> xaxis;
            Vector3<double> yaxis;
            Vector3<double> zaxis;
            if (this->setup_axis(iwater, 0, ito, ith, cel, water[iwater], xaxis, yaxis, zaxis))
            {
                this->wannier_3D(iwater, ito, ith, cel, water[iwater], xaxis, yaxis, zaxis, bond_wan_index, lone_wan_index);
            }
        }

        if(INPUT.func == 2)
        {
            if(fill_shell(cel, water, iwater, first_shell_index, second_shell_index, nwater_first_shell, nwater_second_shell))
            {
                for(int iH=0; iH<2; iH++)
                {
                    Vector3<double> xaxis;
                    Vector3<double> yaxis;
                    Vector3<double> zaxis;
                    if (this->setup_axis(iwater, iH, ito, ith, cel, water[iwater], xaxis, yaxis, zaxis))
                    {
                        this->second_shell_3D(iwater, ito, ith, cel, water, xaxis, yaxis, zaxis, iH, first_shell_index, second_shell_index, nwater_first_shell, nwater_second_shell);
                    }
                }
            }
        }

        if (INPUT.func == 3)// renxi 20220309: calculate the isosurface of first shell neighbors
        {
            if ((INPUT.nacc!=-1 and water[iwater].naccept == INPUT.nacc) or INPUT.nacc == -1)
            {
                if ((INPUT.ndon!=-1 and water[iwater].ndonate == INPUT.ndon) or INPUT.ndon == -1)
                {
                    for(int iH=0; iH<2; iH++)
                    {
                        if (INPUT.ndon != 1 and iH == 1)
                        {
                            continue;
                        }

                        Vector3<double> xaxis;
                        Vector3<double> yaxis;
                        Vector3<double> zaxis;
                        
                        if (this->setup_axis(iwater, iH, ito, ith, cel, water[iwater], xaxis, yaxis, zaxis))
                        {
                            if (INPUT.func_b == 1)
                            {
                                this->first_shell_3D(iwater, ito, ith, cel, water, xaxis, yaxis, zaxis, iH);
                            }
                            else if (INPUT.func_b == 2)
                            {
                                this->second_shell_3D2(iwater, ito, ith, cel, water, xaxis, yaxis, zaxis, iH);
                            }
                        }
                    }
                }
            }
        }


        for(int ish=0; ish<6; ish++)
        {
            delete[] second_shell_index[ish];
        }
        delete[] first_shell_index;
        delete[] second_shell_index;
        delete[] nwater_second_shell;
    }
    delete[] water;

    for (int ia=0; ia<cel.atom[ito].na; ia++)
    {
        delete[] bond_wan_index[ia];
        delete[] lone_wan_index[ia];
    }
    delete[] bond_wan_index;
    delete[] lone_wan_index;

    return;
}

bool Dist2::setup_axis(const int &iwater, const int &iH, const int &ito, const int &ith, const Cell &cel, const Water &water, Vector3<double> &xaxis, Vector3<double> &yaxis, Vector3<double> &zaxis)
{
    // The axises are set as follows: O is the origin. The shorter O-H bond is x axis. The H-O-H plane is the x-y plane. 
    // The y axis points between two O-H bonds.
    // This definition method is biased. The x axis should be chosen symmemetrically.
    //cout << iwater << endl;
    if(water.nH == 2)
    {
        /*int xHindex = -1;
        int Hindex = -1;
        if(water.disH[0] <= water.disH[1])
        {
            xHindex = water.indexH[0];
            Hindex = water.indexH[1];
        }
        else
        {
            xHindex = water.indexH[1];
            Hindex = water.indexH[0];
        }
        Vector3<double> pos_H;
        pos_H.x = cel.atom[ith].pos[xHindex].x;
        pos_H.y = cel.atom[ith].pos[xHindex].y;
        pos_H.z = cel.atom[ith].pos[xHindex].z;
        putback_cell(cel.atom[ito].pos[iwater], pos_H);
        xaxis = pos_H - cel.atom[ito].pos[iwater];
        xaxis = xaxis/xaxis.norm();
        assert(abs(xaxis.norm()-1) <= 1e-12);

        Vector3<double> pos_H1;
        pos_H1.x = cel.atom[ith].pos[Hindex].x;
        pos_H1.y = cel.atom[ith].pos[Hindex].y;
        pos_H1.z = cel.atom[ith].pos[Hindex].z;
        putback_cell(cel.atom[ito].pos[iwater], pos_H1);
        pos_H1 = pos_H1 - cel.atom[ito].pos[iwater];
        pos_H1 = pos_H1/pos_H1.norm();
        double pos_H1_new_x = pos_H1.x*xaxis.x + pos_H1.y*xaxis.y + pos_H1.z*xaxis.z;
        if (pos_H1_new_x >= 0)
        {
            cout << "pos of O " << cel.atom[ito].pos[iwater].x << " " << cel.atom[ito].pos[iwater].y << " " << cel.atom[ito].pos[iwater].z << endl;
            cout << "pos of xH " << pos_H.x << " " << pos_H.y << " " << pos_H.z << endl;
            cout << "pos of H1 " << pos_H1.x << " " << pos_H1.y << " " << pos_H1.z << endl;
            cout << "pos_H1_new_x = " << pos_H1_new_x << endl;
            this->nabnormal++;
            //exit(0);
        }
        yaxis = pos_H1 - pos_H1_new_x*xaxis;
        yaxis = yaxis/yaxis.norm();

        assert(xaxis.x*yaxis.x + xaxis.y*yaxis.y + xaxis.z*yaxis.z <= 1e-12);

        zaxis = xaxis^yaxis;
        assert(xaxis.x*zaxis.x + xaxis.y*zaxis.y + xaxis.z*zaxis.z <= 1e-12);
        assert(yaxis.x*zaxis.x + yaxis.y*zaxis.y + yaxis.z*zaxis.z <= 1e-12);
        return true;*/
        if(INPUT.func == 1)
        {
            Vector3<double> pos_H1;
            Vector3<double> pos_H2;
            pos_H1.x = cel.atom[ith].pos[water.indexH[0]].x;
            pos_H1.y = cel.atom[ith].pos[water.indexH[0]].y;
            pos_H1.z = cel.atom[ith].pos[water.indexH[0]].z;

            pos_H2.x = cel.atom[ith].pos[water.indexH[1]].x;
            pos_H2.y = cel.atom[ith].pos[water.indexH[1]].y;
            pos_H2.z = cel.atom[ith].pos[water.indexH[1]].z;

            putback_cell(cel.atom[ito].pos[iwater], pos_H1);
            putback_cell(cel.atom[ito].pos[iwater], pos_H2);

            pos_H1 = pos_H1 - cel.atom[ito].pos[iwater];
            pos_H2 = pos_H2 - cel.atom[ito].pos[iwater];

            pos_H1 = pos_H1/pos_H1.norm();
            pos_H2 = pos_H2/pos_H2.norm();
            xaxis = (pos_H2 + pos_H1)/(pos_H1+pos_H2).norm();
            double H1x = pos_H1*xaxis;
            yaxis = pos_H1-H1x*xaxis;
            yaxis = yaxis/yaxis.norm();
            assert(abs(xaxis.norm() - 1) <= 1e-12);
            assert(abs(yaxis.norm() - 1) <= 1e-12);
            assert(abs(xaxis*yaxis) <= 1e-12);

            zaxis = xaxis^yaxis;
            assert(abs(zaxis.norm() - 1) <= 1e-12);
            assert(abs(zaxis*yaxis) <= 1e-12);
            assert(abs(xaxis*zaxis) <= 1e-12);

            return true;
        }
        else if (INPUT.func == 2 or INPUT.func == 3)
        {
            int xHindex = -1;
            int Hindex = -1;
            if (iH == 0)
            {
                xHindex = water.indexH[iH];
                Hindex = water.indexH[1];
            }
            else if (iH == 1)
            {
                xHindex = water.indexH[iH];
                Hindex = water.indexH[0];
            }
            else
            {
                cout << "iH input is wrong. iH has to be 0 or 1. " << endl;
                exit(0);
            }
            if (INPUT.ndon == 1) // if ndon=1, select the donated H as x axis
            {
                if (xHindex != water.donateH[0])
                {
                    return false;
                }
                Vector3<double> pos_H;
                pos_H.x = cel.atom[ith].pos[xHindex].x;
                pos_H.y = cel.atom[ith].pos[xHindex].y;
                pos_H.z = cel.atom[ith].pos[xHindex].z;
                putback_cell(cel.atom[ito].pos[iwater], pos_H);
                xaxis = pos_H - cel.atom[ito].pos[iwater];
                xaxis = xaxis/xaxis.norm();
                assert(abs(xaxis.norm()-1) <= 1e-12);

                Vector3<double> pos_H1;
                pos_H1.x = cel.atom[ith].pos[Hindex].x;
                pos_H1.y = cel.atom[ith].pos[Hindex].y;
                pos_H1.z = cel.atom[ith].pos[Hindex].z;
                putback_cell(cel.atom[ito].pos[iwater], pos_H1);
                pos_H1 = pos_H1 - cel.atom[ito].pos[iwater];
                pos_H1 = pos_H1/pos_H1.norm();
                double pos_H1_new_x = pos_H1.x*xaxis.x + pos_H1.y*xaxis.y + pos_H1.z*xaxis.z;
                if (pos_H1_new_x >= 0)
                {
                    cout << "pos of O " << cel.atom[ito].pos[iwater].x << " " << cel.atom[ito].pos[iwater].y << " " << cel.atom[ito].pos[iwater].z << endl;
                    cout << "pos of xH " << pos_H.x << " " << pos_H.y << " " << pos_H.z << endl;
                    cout << "pos of H1 " << pos_H1.x << " " << pos_H1.y << " " << pos_H1.z << endl;
                    cout << "pos_H1_new_x = " << pos_H1_new_x << endl;
                    this->nabnormal++;
                    //exit(0);
                }
                yaxis = pos_H1 - pos_H1_new_x*xaxis;
                yaxis = yaxis/yaxis.norm();

                assert(xaxis.x*yaxis.x + xaxis.y*yaxis.y + xaxis.z*yaxis.z <= 1e-12);

                zaxis = xaxis^yaxis;
                assert(xaxis.x*zaxis.x + xaxis.y*zaxis.y + xaxis.z*zaxis.z <= 1e-12);
                assert(yaxis.x*zaxis.x + yaxis.y*zaxis.y + yaxis.z*zaxis.z <= 1e-12);
                return true;
            }
            else
            {
                Vector3<double> pos_H1 = cel.atom[ith].pos[xHindex];
                Vector3<double> pos_H2 = cel.atom[ith].pos[Hindex];
                Vector3<double> pos_O = cel.atom[ito].pos[iwater];
                putback_cell(pos_O, pos_H1);
                putback_cell(pos_O, pos_H1);
                pos_H1 = (pos_H1 - pos_O)/(pos_H1 - pos_O).norm();
                pos_H2 = (pos_H2 - pos_O)/(pos_H2 - pos_O).norm();
                yaxis = (pos_H1 + pos_H2)/(pos_H1 + pos_H2).norm();
                xaxis = (pos_H1 - pos_H2)/(pos_H1 - pos_H2).norm();
                zaxis = xaxis^yaxis;
                assert(xaxis.x*yaxis.x + xaxis.y*yaxis.y + xaxis.z*yaxis.z <= 1e-12);
                assert(xaxis.x*zaxis.x + xaxis.y*zaxis.y + xaxis.z*zaxis.z <= 1e-12);
                assert(yaxis.x*zaxis.x + yaxis.y*zaxis.y + yaxis.z*zaxis.z <= 1e-12);
                return true;
            }    
        }
    }
    else
    {
        cout << "Number of hydrogen atoms in water is not 2." << endl;
        return false;
    }
}

void Dist2::putback_cell(const Vector3<double> &pos1, Vector3<double> &pos2)
// The function check if pos2 is celldm/2 away from pos1 on every axis. If so, the function will change pos2 according to the periodic boundary condition
// to put it back to pos1's neighborhood.
{
    if(abs(pos1.x-pos2.x) >= INPUT.celldm1/2)
    {
        if (abs(pos1.x - pos2.x -INPUT.celldm1) < abs(pos1.x - pos2.x + INPUT.celldm1) and abs(pos1.x - pos2.x -INPUT.celldm1) <= abs(pos1.x - pos2.x))
        {
            pos2.x += INPUT.celldm1;
        }
        else
        {
            pos2.x -= INPUT.celldm1;
        }
        assert(abs(pos1.x-pos2.x) <= INPUT.celldm1/2);
    }
    if(abs(pos1.y-pos2.y) >= INPUT.celldm2/2)
    {
        if (abs(pos1.y - pos2.y -INPUT.celldm2) < abs(pos1.y - pos2.y + INPUT.celldm2) and abs(pos1.y - pos2.y -INPUT.celldm2) <= abs(pos1.y - pos2.y))
        {
            pos2.y += INPUT.celldm2;
        }
        else
        {
            pos2.y -= INPUT.celldm2;
        }
        assert(abs(pos1.y-pos2.y) <= INPUT.celldm2/2);
    }
    if(abs(pos1.z-pos2.z) >= INPUT.celldm3/2)
    {
        if (abs(pos1.z - pos2.z -INPUT.celldm3) < abs(pos1.z - pos2.z + INPUT.celldm3) and abs(pos1.z - pos2.z -INPUT.celldm3) <= abs(pos1.z - pos2.z))
        {
            pos2.z += INPUT.celldm3;
        }
        else
        {
            pos2.z -= INPUT.celldm3;
        }
        assert(abs(pos1.z-pos2.z) <= INPUT.celldm3/2);
    }
    return;
}

void Dist2::wannier_3D(const int &iwater, const int &ito, const int &ith, const Cell &cel, const Water &water,
 const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis, int** bond_wan_index, int** lone_wan_index)
{
    //cout << "iwater = " << iwater << endl;
    for (int iwan=0; iwan<2; iwan++)
    {
        Vector3<double> wan_pos;
        wan_pos.x = cel.wan_centers[bond_wan_index[iwater][iwan]].x;
        wan_pos.y = cel.wan_centers[bond_wan_index[iwater][iwan]].y;
        wan_pos.z = cel.wan_centers[bond_wan_index[iwater][iwan]].z;

        putback_cell(cel.atom[ito].pos[iwater], wan_pos);
        Vector3<double> new_pos = wan_pos - cel.atom[ito].pos[iwater];
        double x = new_pos*xaxis + INPUT.rcut/2;
        double y = new_pos*yaxis + INPUT.rcut/2;
        double z = new_pos*zaxis + INPUT.rcut/2;
        int which_x = (int) (x/this->dr);
        int which_y = (int) (y/this->dr);
        int which_z = (int) (z/this->dr);
        //cout << "which_x = " << which_x << endl;
        //cout << "which_y = " << which_y << endl;
        //cout << "which_z = " << which_z << endl;
        if (which_x >= INPUT.u1 or which_x < 0)
        {
            cout << "x = " << x << ", which_x = " << which_x << endl;
        }
        if (which_y >= INPUT.u1 or which_y < 0)
        {
            cout << "y = " << y << ", which_y = " << which_y << endl;
        }
        if (which_z >= INPUT.u1 or which_z < 0)
        {
            cout << "z = " << z << ", which_z = " << which_z << endl;
        }
        assert(which_x < INPUT.u1 and which_x >= 0);
        assert(which_y < INPUT.u1 and which_y >= 0);
        assert(which_z < INPUT.u1 and which_z >= 0);
        this->dist3D[which_x][which_y][which_z]++;
    }
    for (int iwan=0; iwan<2; iwan++)
    {
        Vector3<double> wan_pos;
        wan_pos.x = cel.wan_centers[lone_wan_index[iwater][iwan]].x;
        wan_pos.y = cel.wan_centers[lone_wan_index[iwater][iwan]].y;
        wan_pos.z = cel.wan_centers[lone_wan_index[iwater][iwan]].z;

        putback_cell(cel.atom[ito].pos[iwater], wan_pos);
        Vector3<double> new_pos = wan_pos - cel.atom[ito].pos[iwater];
        double x = new_pos*xaxis + INPUT.rcut/2;
        double y = new_pos*yaxis + INPUT.rcut/2;
        double z = new_pos*zaxis + INPUT.rcut/2;
        int which_x = (int) (x/this->dr);
        int which_y = (int) (y/this->dr);
        int which_z = (int) (z/this->dr);
        if (which_x >= INPUT.u1 or which_x < 0)
        {
            cout << "x = " << x << ", which_x = " << which_x << endl;
        }
        if (which_y >= INPUT.u1 or which_y < 0)
        {
            cout << "y = " << y << ", which_y = " << which_y << endl;
        }
        if (which_z >= INPUT.u1 or which_z < 0)
        {
            cout << "z = " << z << ", which_z = " << which_z << endl;
        }
        assert(which_x < INPUT.u1 and which_x >= 0);
        assert(which_y < INPUT.u1 and which_y >= 0);
        assert(which_z < INPUT.u1 and which_z >= 0);
        this->dist3D[which_x][which_y][which_z]++;
    }
    //cout << "iwater = " << iwater << " done." << endl;
    return;
}

bool Dist2::fill_shell(const Cell &cel, const Water* water, const int &iwater, 
int* &first_shell_index, int** &second_shell_index, int &nwater_first_shell, int* &nwater_second_shell)
{
    for(int idon = 0; idon < water[iwater].ndonate; idon++)
    {
        first_shell_index[nwater_first_shell] = water[iwater].donateO[idon];
        nwater_first_shell++;
    }
    for(int iacc = 0; iacc < water[iwater].naccept; iacc++)
    {
        if (nwater_first_shell >= 6)
        {
            cout << "Water has too many HB." << endl;
            return false;
        }
        first_shell_index[nwater_first_shell] = water[iwater].acceptO[iacc];
        nwater_first_shell++;
    }
    for (int iwater2=0; iwater2<nwater_first_shell; iwater2++)
    {
        for(int idon = 0; idon < water[first_shell_index[iwater2]].ndonate; idon++)
        {
            if (water[first_shell_index[iwater2]].donateO[idon] != iwater)
            {
                if(nwater_second_shell[iwater2] >= 6)
                {
                    cout << "Water has too many HB." << endl;
                    return false;
                }
                second_shell_index[iwater2][nwater_second_shell[iwater2]] = water[first_shell_index[iwater2]].donateO[idon];
                nwater_second_shell[iwater2]++;
            }
        }
        for(int iacc = 0; iacc < water[first_shell_index[iwater2]].naccept; iacc++)
        {
            if (water[first_shell_index[iwater2]].acceptO[iacc] != iwater)
            {
                if(nwater_second_shell[iwater2] >= 6)
                {
                    cout << "Water has too many HB." << endl;
                    return false;
                }
                second_shell_index[iwater2][nwater_second_shell[iwater2]] = water[first_shell_index[iwater2]].acceptO[iacc];
                nwater_second_shell[iwater2]++;
            }
        }
    }
    return true;
}

void Dist2::second_shell_3D(const int &iwater, const int &ito, const int &ith, const Cell &cel, Water* &water,
    const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis, const int &iH, 
    int* &first_shell_index, int** &second_shell_index, int &nwater_first_shell, int* &nwater_second_shell)
{
    int H_index = water[iwater].indexH[iH];
    for (int iwater1=0; iwater1<nwater_first_shell; iwater1++)
    {
        for (int iwater2=0; iwater2<nwater_second_shell[iwater1]; iwater2++)
        {
            Vector3<double> pos_O1;
            Vector3<double> pos_O2;
            pos_O1.x = cel.atom[ito].pos[first_shell_index[iwater1]].x;
            pos_O1.y = cel.atom[ito].pos[first_shell_index[iwater1]].y;
            pos_O1.z = cel.atom[ito].pos[first_shell_index[iwater1]].z;

            pos_O2.x = cel.atom[ito].pos[second_shell_index[iwater1][iwater2]].x;
            pos_O2.y = cel.atom[ito].pos[second_shell_index[iwater1][iwater2]].y;
            pos_O2.z = cel.atom[ito].pos[second_shell_index[iwater1][iwater2]].z;

            putback_cell(cel.atom[ito].pos[iwater], pos_O1);
            putback_cell(cel.atom[ito].pos[iwater], pos_O2);

            double HOO_angle = HBs::angle(cel, cel.atom[ith].pos[H_index], cel.atom[ito].pos[iwater], cel.atom[ito].pos[second_shell_index[iwater1][iwater2]]);
            double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[second_shell_index[iwater1][iwater2]], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
            if (HOO_angle >= INPUT.theta_min and HOO_angle <= INPUT.theta and dist >= INPUT.r_min and dist <= INPUT.r_max)
            {
                this->calc_which_HOO_angle(cel, water, ito, ith, iH, iwater, first_shell_index[iwater1], second_shell_index[iwater1][iwater2]);
                double x1 = (pos_O1 - cel.atom[ito].pos[iwater])*xaxis + INPUT.rcut/2;
                double y1 = (pos_O1 - cel.atom[ito].pos[iwater])*yaxis + INPUT.rcut/2;
                double z1 = (pos_O1 - cel.atom[ito].pos[iwater])*zaxis + INPUT.rcut/2;

                double x2 = (pos_O2 - cel.atom[ito].pos[iwater])*xaxis + INPUT.rcut/2;
                double y2 = (pos_O2 - cel.atom[ito].pos[iwater])*yaxis + INPUT.rcut/2;
                double z2 = (pos_O2 - cel.atom[ito].pos[iwater])*zaxis + INPUT.rcut/2;
                
                int which_x1 = (int) (x1/this->dr);
                int which_y1 = (int) (y1/this->dr);
                int which_z1 = (int) (z1/this->dr);

                int which_x2 = (int) (x2/this->dr);
                int which_y2 = (int) (y2/this->dr);
                int which_z2 = (int) (z2/this->dr);

                assert(which_x1 < INPUT.u1 and which_x1 >= 0);
                assert(which_y1 < INPUT.u1 and which_y1 >= 0);
                assert(which_z1 < INPUT.u1 and which_z1 >= 0);

                assert(which_x2 < INPUT.u1 and which_x2 >= 0);
                assert(which_y2 < INPUT.u1 and which_y2 >= 0);
                assert(which_z2 < INPUT.u1 and which_z2 >= 0);

                this->dist3D[which_x1][which_y1][which_z1]++;
                this->dist3D2[which_x2][which_y2][which_z2]++;
            }
        }
    }
    return;
}

void Dist2::calc_which_HOO_angle(const Cell &cel, Water* &water, const int &ito, const int &ith, const int &iH, const int &water_index, int &water_index1, int &water_index2)
{
    int iH2 = (int) abs(1-iH);
    double angle1 = HBs::angle(cel, cel.atom[ito].pos[water_index1], cel.atom[ito].pos[water_index], cel.atom[ith].pos[water[water_index].indexH[iH]]);
    double angle2 = HBs::angle(cel, cel.atom[ito].pos[water_index1], cel.atom[ito].pos[water_index], cel.atom[ith].pos[water[water_index].indexH[iH2]]);
    bool thisHdonate = false;
    bool otherHdonate = false;
    bool seconddonate = false;
    bool secondaccept = false;
    if(angle1 < INPUT.acut_hoo)
    {
        thisHdonate = true;
    }
    else if (angle2 < INPUT.acut_hoo)
    {
        otherHdonate = true;
    }

    for (int idon=0; idon<water[water_index1].ndonate; idon++)
    {
        if (water[water_index1].donateO[idon] == water_index2)
        {
            seconddonate = true;
            break;
        }
    }
    if (!seconddonate)
    {
        for (int iacc=0; iacc<water[water_index1].naccept; iacc++)
        {
            secondaccept = true;
            break;
        }
    }
    if (secondaccept == false and seconddonate == false)
    {
        cout << "Second and first shell is not connected. Check." << endl;
        exit(0);
    }
    if (thisHdonate and seconddonate)
    {
        this->which_HOO_angle[0]++;
    }
    else if (otherHdonate and seconddonate)
    {
        this->which_HOO_angle[1]++;
    }
    else if (seconddonate)
    {
        this->which_HOO_angle[2]++;
    }
    else if (thisHdonate and secondaccept)
    {
        this->which_HOO_angle[3]++;
    }
    else if (otherHdonate and secondaccept)
    {
        this->which_HOO_angle[4]++;
    }
    else
    {
        this->which_HOO_angle[5]++;
    }
    return;
}

void Dist2::first_shell_3D(const int &iwater, const int &ito, const int &ith, const Cell &cel, Water* &water, 
    const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis, const int &iH)
{
    for (int idon=0; idon<water[iwater].ndonate; idon++)
    {
        int Oindex = water[iwater].donateO[idon];
        this->addup_3D(cel, ito, iwater, Oindex, xaxis, yaxis, zaxis);
    }
    for (int iacc=0; iacc<water[iwater].naccept; iacc++)
    {
        int Oindex = water[iwater].acceptO[iacc];
        this->addup_3D(cel, ito, iwater, Oindex, xaxis, yaxis, zaxis);
    }
    return;
}

void Dist2::second_shell_3D2(const int &iwater, const int &ito, const int &ith, const Cell &cel, Water* &water, 
    const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis, const int &iH)
{
    int* first_shell_index = new int[10];
    int nwater_in_first_shell = 0;

    for (int idon=0; idon<water[iwater].ndonate; idon++)
    {
        first_shell_index[idon] = water[iwater].donateO[idon];
    }
    for (int iacc=0; iacc<water[iwater].naccept; iacc++)
    {
        first_shell_index[iacc+water[iwater].ndonate] = water[iwater].acceptO[iacc];
    }
    nwater_in_first_shell = water[iwater].ndonate + water[iwater].naccept;

    for (int iwater2=0; iwater2<cel.atom[ito].na; iwater2++)
    {
        if (iwater2 == iwater) return;
        for (int iwater3=0; iwater3<nwater_in_first_shell; iwater3++)
        {
            if(first_shell_index[iwater3] == iwater2) return;
        }
        if (distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3) <= INPUT.rcut1)
        {
            this->addup_3D(cel, ito, iwater, iwater2, xaxis, yaxis, zaxis);
        }
    }
    return;
}

void Dist2::addup_3D(const Cell &cel, const int &ito, const int &iwater, const int &Oindex, const Vector3<double> &xaxis, const Vector3<double> &yaxis, const Vector3<double> &zaxis)
{
    Vector3<double> pos_O = cel.atom[ito].pos[Oindex];
    putback_cell(cel.atom[ito].pos[iwater], pos_O);
    Vector3<double> delta_pos = pos_O - cel.atom[ito].pos[iwater];
    double new_x = delta_pos*xaxis + INPUT.rcut/2;
    double new_y = delta_pos*yaxis + INPUT.rcut/2;
    double new_z = delta_pos*zaxis + INPUT.rcut/2;
    if (INPUT.func_c == 2 and new_z >= INPUT.rcut/2)
    {
        return;
    }
    if (new_x < INPUT.rcut and new_y < INPUT.rcut and new_z < INPUT.rcut and new_x > 0 and new_y > 0 and new_z > 0)
    {
        int whichx = (int) (new_x/this->dr);
        int whichy = (int) (new_y/this->dr);
        int whichz = (int) (new_z/this->dr);
        this->dist3D[whichx][whichy][whichz]++;
    }
    return;
}