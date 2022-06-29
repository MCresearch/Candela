#include "cellFile.h"
#include "input.h"
#include "incremental_pdf2.h"
#include "math.h"
#include "HBs.h"
#include "gfun.h"
#include "pdf_added.h"

incrementalPDF2::incrementalPDF2(){}
incrementalPDF2::~incrementalPDF2(){}
void incrementalPDF2::Routine()
{
	TITLE("incrementalPDF2", "Routine");
	cout << "Calculate the radial distribution functions g(r)." << endl;

	assert(INPUT.dr>0.0);

	assert(INPUT.rcut>0.0);

    assert(INPUT.nshell > 0);

    this->allocate();

	int count_geometry_number = 0;
	for(int igeo=INPUT.geo_1; igeo<=INPUT.geo_2; ++igeo) 
	{
		cout << "igeo=" << igeo << endl;
		CellFile cel;

		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;
		cout << "Succeeded" << endl;
		stringstream ss; ss << igeo;
		cel.file_name = ss.str();
		// cel : input geometry file
		if( !CellFile::ReadGeometry( cel ) ) continue;

		if(cel.read_and_used==false) 
		{
			cel.clean();//qianrui add
			continue;
		}

        this->calc(cel);

        cel.clean();
        count_geometry_number++;
        cout << "count_geometry_number " << count_geometry_number << endl;

        if (INPUT.func == 2 and this->ngeometry >= INPUT.ntry)
        {
            break;
        }

    }//igeo 

    double rho = INPUT.natom1/INPUT.celldm1/INPUT.celldm2/INPUT.celldm3;
    for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
    {
        for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
        {
            this->gr[ishell][ir] /= double(count_geometry_number)*rho*4*PI*pow((ir+0.5), 2)*pow(INPUT.dr, 3)*INPUT.natom1;
        }
    }

    this->out();
    return;
}

void incrementalPDF2::allocate()
{
    int expnshell = pow(2, INPUT.nshell);
    this->gr = new double*[expnshell];
    // for nshell = 2, 0: D->D; 1: D->A; 2: A->D; 3: A->A;
    // for nshell = 3, 0: D->D->D; 1: D->D->A; 2: D->A->D; 3: D->A->A; 4: A->D->D; 5: A->D->A; 6: A->A->D; 7: A->A->A.
    this->calculated = new bool*[expnshell];
    for (int ishell=0; ishell<expnshell; ishell++)
    {
        this->gr[ishell] = new double[int(INPUT.rcut/INPUT.dr)+1];
        this->calculated[ishell] = new bool[INPUT.natom1];
        for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
        {
            gr[ishell][ir] = 0;
        }
        for (int ia=0; ia<INPUT.natom1; ia++)
        {
            calculated[ishell][ia] = false;
        }
    }
    if (INPUT.nshell == 3)
    {
        this->angle_planar_total = new double[int(INPUT.theta/INPUT.dtheta)+1];
        this->angle_planar = new double*[expnshell];
        this->angle_planar_OOO = new double**[expnshell];
        this->angle_planar_OOO_total = new double*[int(INPUT.theta/INPUT.dtheta)+1];
        for (int ishell=0; ishell<expnshell; ishell++)
        {
            this->angle_planar[ishell] = new double[int(INPUT.theta/INPUT.dtheta)+1];
            this->angle_planar_OOO[ishell] = new double*[int(INPUT.theta/INPUT.dtheta)+1];
            for (int ir=0; ir<int(INPUT.theta/INPUT.dtheta)+1; ir++)
            {
                this->angle_planar[ishell][ir] = 0;
                this->angle_planar_OOO[ishell][ir] = new double[2*int(INPUT.theta/INPUT.dtheta)+1];
                for (int ir2=0; ir2<2*int(INPUT.theta/INPUT.dtheta)+1; ir2++)
                {
                    this->angle_planar_OOO[ishell][ir][ir2] = 0;
                }
            }
        }
        for (int ir=0; ir<int(INPUT.theta/INPUT.dtheta)+1; ir++)
        {
            this->angle_planar_total[ir] = 0;
            this->angle_planar_OOO_total[ir] = new double[2*int(INPUT.theta/INPUT.dtheta)+1];
            for (int ir2=0; ir2<2*int(INPUT.theta/INPUT.dtheta)+1; ir2++)
            {
                this->angle_planar_OOO_total[ir][ir2] = 0;
            }
        }

        this->distance_12_34_14_total = new double*[int(INPUT.rcut/INPUT.dr)+1];
        this->distance_12_34_14 = new double**[expnshell];
        for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
        {
            this->distance_12_34_14_total[ir] = new double[2*int(INPUT.rcut/INPUT.dr)+1];
            for (int ir2=0; ir2<2*int(INPUT.rcut/INPUT.dr)+1; ir2++)
            {
                this->distance_12_34_14_total[ir][ir2] = 0;
            }
        }
        for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
        {
            this->distance_12_34_14[ishell] = new double*[int(INPUT.rcut/INPUT.dr)+1];
            for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
            {
                this->distance_12_34_14[ishell][ir] = new double[2*int(INPUT.rcut/INPUT.dr)+1];
                for (int ir2=0; ir2<2*int(INPUT.rcut/INPUT.dr)+1; ir2++)
                {
                    this->distance_12_34_14[ishell][ir][ir2] = 0;
                }
            }
        }

        this->angle_OOO_2_total = new double*[int(INPUT.theta/INPUT.dtheta)+1];
        this->angle_OOO_2 = new double**[expnshell];
        for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta)+1; ia++)
        {
            this->angle_OOO_2_total[ia] = new double[int(INPUT.theta/INPUT.dtheta)+1];
            for (int ia2=0; ia2<int(INPUT.theta/INPUT.dtheta)+1; ia2++)
            {
                this->angle_OOO_2_total[ia][ia2] = 0;
            }
        }
        for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
        {
            this->angle_OOO_2[ishell] = new double*[int(INPUT.theta/INPUT.dtheta)+1];
            for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta)+1; ia++)
            {
                this->angle_OOO_2[ishell][ia] = new double[int(INPUT.theta/INPUT.dtheta)+1];
                for (int ia2=0; ia2<int(INPUT.theta/INPUT.dtheta)+1; ia2++)
                {
                    this->angle_OOO_2[ishell][ia][ia2] = 0;
                }
            }
        }
        if (INPUT.func == 2)
        {
            this->ofs_xyz1.open("third_shell_geometry_1.xyz");
            this->ofs_xyz2.open("third_shell_geometry_2.xyz");
            this->ofs_xyz3.open("third_shell_geometry_3.xyz");
            this->ofs_xyz4.open("third_shell_geometry_4.xyz");
            this->ngeometry = 0;
        }
    }
    return;
}

void incrementalPDF2::set_false(Cell &cel, Water* water, int &iwater)
{
    for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
    {
        for (int ia=0; ia<INPUT.natom1; ia++)
        {
            this->calculated[ishell][ia] = false;
        }
        this->calculated[ishell][iwater] = true;
    }
    
    if (INPUT.nshell >= 2)
    {
        for (int idon=0; idon<water[iwater].ndonate; idon++)
        {
            int iwater2 = water[iwater].donateO[idon];
            for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
            {
                this->calculated[ishell][iwater2] = true;
            }
            if (INPUT.nshell >= 3)
            {
                for (int idon2=0; idon2<water[iwater2].ndonate; idon2++)
                {
                    int iwater3=water[iwater2].donateO[idon2];
                    for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
                    {
                        this->calculated[ishell][iwater3] = true;
                    }
                }
                for (int iacc2=0; iacc2<water[iwater2].naccept; iacc2++)
                {
                    int iwater3=water[iwater2].acceptO[iacc2];
                    for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
                    {
                        this->calculated[ishell][iwater3] = true;
                    }
                }
            }
        }
        for (int iacc=0; iacc<water[iwater].naccept; iacc++)
        {
            int iwater2 = water[iwater].acceptO[iacc];
            for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
            {
                this->calculated[ishell][iwater2] = true;
            }
            if (INPUT.nshell >= 3)
            {
                for (int idon2=0; idon2<water[iwater2].ndonate; idon2++)
                {
                    int iwater3=water[iwater2].donateO[idon2];
                    for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
                    {
                        this->calculated[ishell][iwater3] = true;
                    }
                }
                for (int iacc2=0; iacc2<water[iwater2].naccept; iacc2++)
                {
                    int iwater3=water[iwater2].acceptO[iacc2];
                    for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
                    {
                        this->calculated[ishell][iwater3] = true;
                    }
                }
            }
        }
    }
    return;
}

void incrementalPDF2::calc(Cell &cel)
{
	int ito = -1;
	int ith = -1;
	int itc = -1;
	int itcl = -1;
	int itna = -1;
	int it_ele1 = -1;
	int it_ele2 = -1;
	for(int it=0;it <INPUT.ntype; ++it)
	{
		if(cel.atom[it].id=="O") ito=it;
		else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;    
    }

    Water *water;
	assert(INPUT.system == "water");
	water = new Water[cel.atom[ito].na];
	Water::nions = 0;
	HBs::setup_water(cel, water);
    if (INPUT.nshell == 2)
    {
        for (int iwater=0; iwater<cel.atom[ito].na; iwater++)
        {
            this->set_false(cel, water, iwater);
            for (int idon=0; idon<water[iwater].ndonate; idon++)
            {
                int iwater2 = water[iwater].donateO[idon];
                for (int idon2=0; idon2<water[iwater2].ndonate; idon2++)
                {
                    int iwater3 = water[iwater2].donateO[idon2];
                    if (iwater3 != iwater and this->calculated[0][iwater3] == false)
                    {
                        double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater3], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                        if (dist < INPUT.rcut)
                        {
                            int which = int(dist/INPUT.dr);
                            this->gr[0][which]++;
                            this->calculated[0][iwater3] = true;
                        }
                    }
                }

                for (int iacc2=0; iacc2<water[iwater2].naccept; iacc2++)
                {
                    int iwater3 = water[iwater2].acceptO[iacc2];
                    if (iwater3 != iwater and this->calculated[1][iwater3] == false)
                    {
                        double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater3], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                        if (dist < INPUT.rcut)
                        {
                            int which = int(dist/INPUT.dr);
                            this->gr[1][which]++;
                            this->calculated[1][iwater3] = true;
                        }
                    }
                }
            }

            for (int iacc=0; iacc<water[iwater].naccept; iacc++)
            {
                int iwater2 = water[iwater].acceptO[iacc];
                for (int idon2=0; idon2<water[iwater2].ndonate; idon2++)
                {
                    int iwater3 = water[iwater2].donateO[idon2];
                    if (iwater3 != iwater and this->calculated[2][iwater3] == false)
                    {
                        double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater3], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                        if (dist < INPUT.rcut)
                        {
                            int which = int(dist/INPUT.dr);
                            this->gr[2][which]++;
                            this->calculated[2][iwater3] = true;
                        }
                    }
                }

                for (int iacc2=0; iacc2<water[iwater2].naccept; iacc2++)
                {
                    int iwater3 = water[iwater2].acceptO[iacc2];
                    if (iwater3 != iwater and this->calculated[3][iwater3] == false)
                    {
                        double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater3], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                        if (dist < INPUT.rcut)
                        {
                            int which = int(dist/INPUT.dr);
                            this->gr[3][which]++;
                            this->calculated[3][iwater3] = true;
                        }
                    }
                }
            }
        }
    }
    else if (INPUT.nshell == 3)
    {
        for (int iwater=0; iwater<cel.atom[ito].na; iwater++)
        {
            this->set_false(cel, water, iwater);
            for (int idon=0; idon<water[iwater].ndonate; idon++)
            {
                int iwater2 = water[iwater].donateO[idon];
                for (int idon2=0; idon2<water[iwater2].ndonate; idon2++)
                {
                    int iwater3 = water[iwater2].donateO[idon2];
                    if (iwater3 != iwater)
                    {
                        for (int idon3=0; idon3<water[iwater3].ndonate; idon3++)
                        {
                            int iwater4 = water[iwater3].donateO[idon3];
                            if (iwater4 != iwater and iwater4 != iwater2 and calculated[0][iwater4] == false)
                            {
                                double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                if (dist < INPUT.rcut)
                                {
                                    /*if (dist < 3.1) 
                                    {
                                        cout << "DDD " << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << " " << dist << endl;
                                    }*/
                                    double angle = this->calc_angle_planar(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith);
                                    double angle_ooo1 = HBs::angle(cel, cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3]);
                                    double angle_ooo2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4]);
                                    this->angle_planar_OOO_total[int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_OOO[0][int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_total[int(angle/INPUT.dtheta)]++;
                                    this->angle_planar[0][int(angle/INPUT.dtheta)]++;
                                    int which = int(dist/INPUT.dr);
                                    this->gr[0][which]++;
                                    calculated[0][iwater4] = true;
                                    this->angle_OOO_2_total[int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    this->angle_OOO_2[0][int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    double dist1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    double dist2 = distance(cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    this->distance_12_34_14[0][int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    this->distance_12_34_14_total[int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    if (INPUT.func == 2)
                                    {
                                        this->output_geometry(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith, 0);
                                    }
                                }
                            }
                        }
                        for (int iacc3=0; iacc3<water[iwater3].naccept; iacc3++)
                        {
                            int iwater4 = water[iwater3].acceptO[iacc3];
                            if (iwater4 != iwater and iwater4 != iwater2 and calculated[1][iwater4] == false)
                            {
                                double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                if (dist < INPUT.rcut)
                                {
                                    /*if (dist < 3.1) 
                                    {
                                        cout << "DDA " << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << " " << dist << endl;
                                    }*/
                                    double angle = this->calc_angle_planar(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith);
                                    double angle_ooo1 = HBs::angle(cel, cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3]);
                                    double angle_ooo2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4]);
                                    this->angle_planar_OOO_total[int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_OOO[1][int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_total[int(angle/INPUT.dtheta)]++;
                                    this->angle_planar[1][int(angle/INPUT.dtheta)]++;
                                    int which = int(dist/INPUT.dr);
                                    this->gr[1][which]++;
                                    calculated[1][iwater4] == true;
                                    this->angle_OOO_2_total[int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    this->angle_OOO_2[1][int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    double dist1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    double dist2 = distance(cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    this->distance_12_34_14[1][int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    this->distance_12_34_14_total[int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    if (INPUT.func == 2)
                                    {
                                        this->output_geometry(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith, 1);
                                    }
                                }
                            }
                        }
                    }
                }
                for (int iacc2=0; iacc2<water[iwater2].naccept; iacc2++)
                {
                    int iwater3 = water[iwater2].acceptO[iacc2];
                    if (iwater3 != iwater)
                    {
                        for (int idon3=0; idon3<water[iwater3].ndonate; idon3++)
                        {
                            int iwater4 = water[iwater3].donateO[idon3];
                            if (iwater4 != iwater and iwater4 != iwater2 and calculated[2][iwater4] == false)
                            {
                                double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                if (dist < INPUT.rcut)
                                {
                                    /*if (dist < 3.1) 
                                    {
                                        cout << "DAD " << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << " " << dist << endl;
                                    }*/
                                    double angle = this->calc_angle_planar(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith);
                                    double angle_ooo1 = HBs::angle(cel, cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3]);
                                    double angle_ooo2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4]);
                                    this->angle_planar_OOO_total[int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_OOO[2][int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_total[int(angle/INPUT.dtheta)]++;
                                    this->angle_planar[2][int(angle/INPUT.dtheta)]++;
                                    int which = int(dist/INPUT.dr);
                                    this->gr[2][which]++;
                                    calculated[2][iwater4] == true;
                                    this->angle_OOO_2_total[int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    this->angle_OOO_2[2][int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    double dist1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    double dist2 = distance(cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    this->distance_12_34_14[2][int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    this->distance_12_34_14_total[int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    if (INPUT.func == 2)
                                    {
                                        this->output_geometry(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith, 2);
                                    }
                                }
                            }
                        }
                        for (int iacc3=0; iacc3<water[iwater3].naccept; iacc3++)
                        {
                            int iwater4 = water[iwater3].acceptO[iacc3];
                            if (iwater4 != iwater and iwater4 != iwater2 and calculated[3][iwater4] == false)
                            {
                                double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                if (dist < INPUT.rcut)
                                {
                                    /*if (dist < 3.1) 
                                    {
                                        cout << "DAA " << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << " " << dist << endl;
                                    }*/
                                    double angle = this->calc_angle_planar(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith);
                                    double angle_ooo1 = HBs::angle(cel, cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3]);
                                    double angle_ooo2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4]);
                                    this->angle_planar_OOO_total[int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_OOO[3][int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_total[int(angle/INPUT.dtheta)]++;
                                    this->angle_planar[3][int(angle/INPUT.dtheta)]++;
                                    int which = int(dist/INPUT.dr);
                                    this->gr[3][which]++;
                                    calculated[3][iwater4] == true;
                                    this->angle_OOO_2_total[int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    this->angle_OOO_2[3][int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    double dist1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    double dist2 = distance(cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    this->distance_12_34_14[3][int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    this->distance_12_34_14_total[int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    if (INPUT.func == 2)
                                    {
                                        this->output_geometry(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith, 3);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for (int iacc=0; iacc<water[iwater].naccept; iacc++)
            {
                int iwater2 = water[iwater].acceptO[iacc];
                for (int idon2=0; idon2<water[iwater2].ndonate; idon2++)
                {
                    int iwater3 = water[iwater2].donateO[idon2];
                    if (iwater3 != iwater)
                    {
                        for (int idon3=0; idon3<water[iwater3].ndonate; idon3++)
                        {
                            int iwater4 = water[iwater3].donateO[idon3];
                            if (iwater4 != iwater and iwater4 != iwater2 and calculated[4][iwater4] == false)
                            {
                                double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                if (dist < INPUT.rcut)
                                {
                                    /*if (dist < 3.1) 
                                    {
                                        cout << "ADD " << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << " " << dist << endl;
                                    }*/
                                    double angle = this->calc_angle_planar(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith);
                                    double angle_ooo1 = HBs::angle(cel, cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3]);
                                    double angle_ooo2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4]);
                                    this->angle_planar_OOO_total[int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_OOO[4][int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_total[int(angle/INPUT.dtheta)]++;
                                    this->angle_planar[4][int(angle/INPUT.dtheta)]++;
                                    int which = int(dist/INPUT.dr);
                                    this->gr[4][which]++;
                                    calculated[4][iwater4] == true;
                                    this->angle_OOO_2_total[int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    this->angle_OOO_2[4][int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    double dist1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    double dist2 = distance(cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    this->distance_12_34_14[4][int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    this->distance_12_34_14_total[int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    if (INPUT.func == 2)
                                    {
                                        this->output_geometry(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith, 4);
                                    }
                                }
                            }
                        }
                        for (int iacc3=0; iacc3<water[iwater3].naccept; iacc3++)
                        {
                            int iwater4 = water[iwater3].acceptO[iacc3];
                            if (iwater4 != iwater and iwater4 != iwater2 and calculated[5][iwater4] == false)
                            {
                                double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                if (dist < INPUT.rcut)
                                {
                                    /*if (dist < 3.1) 
                                    {
                                        cout << "ADA " << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << " " << dist << endl;
                                    }*/
                                    double angle = this->calc_angle_planar(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith);
                                    double angle_ooo1 = HBs::angle(cel, cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3]);
                                    double angle_ooo2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4]);
                                    this->angle_planar_OOO_total[int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_OOO[5][int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_total[int(angle/INPUT.dtheta)]++;
                                    this->angle_planar[5][int(angle/INPUT.dtheta)]++;
                                    int which = int(dist/INPUT.dr);
                                    this->gr[5][which]++;
                                    calculated[5][iwater4] == true;
                                    this->angle_OOO_2_total[int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    this->angle_OOO_2[5][int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    double dist1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    double dist2 = distance(cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    this->distance_12_34_14[5][int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    this->distance_12_34_14_total[int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    if (INPUT.func == 2)
                                    {
                                        this->output_geometry(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith, 5);
                                    }
                                }
                            }
                        }
                    }
                }
                for (int iacc2=0; iacc2<water[iwater2].naccept; iacc2++)
                {
                    int iwater3 = water[iwater2].acceptO[iacc2];
                    if (iwater3 != iwater)
                    {
                        for (int idon3=0; idon3<water[iwater3].ndonate; idon3++)
                        {
                            int iwater4 = water[iwater3].donateO[idon3];
                            if (iwater4 != iwater and iwater4 != iwater2 and calculated[6][iwater4] == false)
                            {
                                double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                if (dist < INPUT.rcut)
                                {
                                    /*if (dist < 3.1) 
                                    {
                                        cout << "AAD " << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << " " << dist << endl;
                                    }*/
                                    double angle = this->calc_angle_planar(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith);
                                    double angle_ooo1 = HBs::angle(cel, cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3]);
                                    double angle_ooo2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4]);
                                    this->angle_planar_OOO_total[int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_OOO[6][int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_total[int(angle/INPUT.dtheta)]++;
                                    this->angle_planar[6][int(angle/INPUT.dtheta)]++;
                                    int which = int(dist/INPUT.dr);
                                    this->gr[6][which]++;
                                    calculated[6][iwater4] == true;
                                    this->angle_OOO_2_total[int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    this->angle_OOO_2[6][int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    double dist1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    double dist2 = distance(cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    this->distance_12_34_14[6][int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    this->distance_12_34_14_total[int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    if (INPUT.func == 2)
                                    {
                                        this->output_geometry(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith, 6);
                                    }
                                }
                            }
                        }
                        for (int iacc3=0; iacc3<water[iwater3].naccept; iacc3++)
                        {
                            int iwater4 = water[iwater3].acceptO[iacc3];
                            if (iwater4 != iwater and iwater4 != iwater2 and calculated[7][iwater4] == false)
                            {
                                double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                if (dist < INPUT.rcut)
                                {
                                    /*if (dist < 3.1) 
                                    {
                                        cout << "AAA " << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << " " << dist << endl;
                                    }*/
                                    double angle = this->calc_angle_planar(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith);
                                    double angle_ooo1 = HBs::angle(cel, cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3]);
                                    double angle_ooo2 = HBs::angle(cel, cel.atom[ito].pos[iwater2], cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4]);
                                    this->angle_planar_OOO_total[int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_OOO[7][int(angle/INPUT.dtheta)][int((angle_ooo1+angle_ooo2)/INPUT.dtheta)]++;
                                    this->angle_planar_total[int(angle/INPUT.dtheta)]++;
                                    this->angle_planar[7][int(angle/INPUT.dtheta)]++;
                                    int which = int(dist/INPUT.dr);
                                    this->gr[7][which]++;
                                    calculated[7][iwater4] == true;
                                    this->angle_OOO_2_total[int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    this->angle_OOO_2[7][int(angle_ooo1/INPUT.dtheta)][int(angle_ooo2/INPUT.dtheta)]++;
                                    double dist1 = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    double dist2 = distance(cel.atom[ito].pos[iwater3], cel.atom[ito].pos[iwater4], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                                    this->distance_12_34_14[7][int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    this->distance_12_34_14_total[int(dist/INPUT.dr)][int((dist1+dist2)/INPUT.dr)]++;
                                    if (INPUT.func == 2)
                                    {
                                        this->output_geometry(cel, water, iwater, iwater2, iwater3, iwater4, ito, ith, 7);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    delete[] water;
    return;
}

void incrementalPDF2::out()
{
    ofstream ofs("pdf.txt");
    for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
    {
        ofs << (ir+0.5)*INPUT.dr << " ";
        for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
        {
            ofs << gr[ishell][ir] << " ";
        }
        ofs << endl;
    }
    ofs.close();
    if (INPUT.nshell == 3)
    {
        ofstream ofs_angle("planar_angle.txt");
        double sum = 0;
        for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta); ia++)
        {
            sum += this->angle_planar_total[ia];
        }
        sum *= INPUT.dtheta;
        for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta); ia++)
        {
            ofs_angle << (ia+0.5)*INPUT.dtheta << " ";
            for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
            {
                ofs_angle << this->angle_planar[ishell][ia]/sum << " ";
            }
            ofs_angle << this->angle_planar_total[ia]/sum << endl;
        }
        ofs_angle.close();

        ofstream ofs_angle_2D("angle_planar_OOO.txt");
        sum = 0;
        for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta); ia++)
        {
            for (int ia2=0; ia2<int(INPUT.theta*2/INPUT.dtheta); ia2++)
            {
                sum += this->angle_planar_OOO_total[ia][ia2];
            }
        }
        sum *= INPUT.dtheta*INPUT.dtheta;
        ofs_angle_2D << "  ";
        for (int ia=0; ia<int(INPUT.theta*2/INPUT.dtheta); ia++)
        {
            ofs_angle_2D << ia*INPUT.dtheta << " ";
        }
        ofs_angle_2D << endl;
        for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta); ia++)
        {
            ofs_angle_2D << ia*INPUT.dtheta << " ";
            for (int ia2=0; ia2<int(2*INPUT.theta/INPUT.dtheta); ia2++)
            {
                ofs_angle_2D << this->angle_planar_OOO_total[ia][ia2]/sum << " ";
            }
            ofs_angle_2D << endl;
        }
        ofs_angle_2D.close();
    
        for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
        {
            ofstream ofs_angle_2D_seperate("angle_planar_OOO_" + to_string(ishell+1) + ".txt");
            sum = 0;
            for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta); ia++)
            {
                for (int ia2=0; ia2<int(INPUT.theta*2/INPUT.dtheta); ia2++)
                {
                    sum += this->angle_planar_OOO[ishell][ia][ia2];
                }
            }
            sum *= INPUT.dtheta*INPUT.dtheta;
            ofs_angle_2D_seperate << "  ";
            for (int ia=0; ia<int(INPUT.theta*2/INPUT.dtheta); ia++)
            {
                ofs_angle_2D_seperate << ia*INPUT.dtheta << " ";
            }
            ofs_angle_2D_seperate << endl;
            for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta); ia++)
            {
                ofs_angle_2D_seperate << ia*INPUT.dtheta << " ";
                for (int ia2=0; ia2<int(2*INPUT.theta/INPUT.dtheta); ia2++)
                {
                    ofs_angle_2D_seperate << this->angle_planar_OOO[ishell][ia][ia2]/sum << " ";
                }
                ofs_angle_2D_seperate << endl;
            }
            ofs_angle_2D_seperate.close();
        }

        ofstream ofs_OOO2_total("OOO2.txt");
        sum = 0;
        for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta)+1; ia++)
        {
            for (int ia2=0; ia2<int(INPUT.theta/INPUT.dtheta)+1; ia2++)
            {
                sum += this->angle_OOO_2_total[ia][ia2];
            }
        }
        sum *= INPUT.dtheta*INPUT.dtheta;
        ofs_OOO2_total << "  ";
        for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta)+1; ia++)
        {
            ofs_OOO2_total << (ia+0.5)*INPUT.dtheta << " ";
        }
        ofs_OOO2_total << endl;
        for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta)+1; ia++)
        {
            ofs_OOO2_total << (ia+0.5)*INPUT.dtheta << " ";
            for (int ia2=0; ia2<int(INPUT.theta/INPUT.dtheta)+1; ia2++)
            {
                ofs_OOO2_total << this->angle_OOO_2_total[ia][ia2]/sum << " ";
            }
            ofs_OOO2_total << endl;
        }
        ofs_OOO2_total.close();
        for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
        {
            ofstream ofs_OOO2("OOO2_" + to_string(ishell) + ".txt");
            sum = 0;
            for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta)+1; ia++)
            {
                for (int ia2=0; ia2<int(INPUT.theta/INPUT.dtheta)+1; ia2++)
                {
                    sum += this->angle_OOO_2[ishell][ia][ia2];
                }
            }
            sum *= INPUT.dtheta*INPUT.dtheta;
            ofs_OOO2 << "  ";
            for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta)+1; ia++)
            {
                ofs_OOO2 << (ia+0.5)*INPUT.dtheta << " ";
            }
            ofs_OOO2 << endl;
            for (int ia=0; ia<int(INPUT.theta/INPUT.dtheta)+1; ia++)
            {
                ofs_OOO2 << (ia+0.5)*INPUT.dtheta << " ";
                for (int ia2=0; ia2<int(INPUT.theta/INPUT.dtheta)+1; ia2++)
                {
                    ofs_OOO2 << this->angle_OOO_2[ishell][ia][ia2]/sum << " ";
                }
                ofs_OOO2 << endl;
            }
            ofs_OOO2.close();
        }

        ofstream ofs_distance("distance_12_34_14_total.txt");
        sum = 0;
        for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
        {
            for (int ir2=0; ir2<int(INPUT.rcut/INPUT.dr)*2+1; ir2++)
            {
                sum += this->distance_12_34_14_total[ir][ir2];
            }
        }
        sum *= INPUT.dr*INPUT.dr;
        ofs_distance << "  ";
        for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)*2+1; ir++)
        {
            ofs_distance << (ir+0.5)*INPUT.dr << " ";
        }
        ofs_distance << endl;
        for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
        {
            ofs_distance << (ir+0.5)*INPUT.dr << " ";
            for (int ir2=0; ir2<int(INPUT.rcut/INPUT.dr)*2+1; ir2++)
            {
                ofs_distance << this->distance_12_34_14_total[ir][ir2]/sum << " ";
            }
            ofs_distance << endl;
        }
        ofs_distance.close();

        for (int ishell=0; ishell<pow(2, INPUT.nshell); ishell++)
        {
            ofstream ofs_distance_sep("distance_12_34_14_" + to_string(ishell) + ".txt");
            sum = 0;
            for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
            {
                for (int ir2=0; ir2<int(INPUT.rcut/INPUT.dr)*2+1; ir2++)
                {
                    sum += this->distance_12_34_14[ishell][ir][ir2];
                }
            }
            sum *= INPUT.dr*INPUT.dr;
            ofs_distance_sep << "  ";
            for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)*2+1; ir++)
            {
                ofs_distance_sep << (ir+0.5)*INPUT.dr << " ";
            }
            ofs_distance_sep << endl;
            for (int ir=0; ir<int(INPUT.rcut/INPUT.dr)+1; ir++)
            {
                ofs_distance_sep << (ir+0.5)*INPUT.dr << " ";
                for (int ir2=0; ir2<int(INPUT.rcut/INPUT.dr)*2+1; ir2++)
                {
                    ofs_distance_sep << this->distance_12_34_14[ishell][ir][ir2]/sum << " ";
                }
                ofs_distance_sep << endl;
            }
            ofs_distance_sep.close();
        }
        if (INPUT.func == 2)
        {
            this->ofs_xyz1.close();
            this->ofs_xyz2.close();
            this->ofs_xyz3.close();
            this->ofs_xyz4.close();
        }
    }// if nshell==3
    return;
}

double incrementalPDF2::calc_angle_planar(Cell &cel, Water* water, int &iwater, int &iwater2, int &iwater3, int &iwater4, int &ito, int &ith)
{
    Vector3<double> pos1;
    Vector3<double> pos2;
    Vector3<double> pos3;
    Vector3<double> pos4;
    pos1.x = cel.atom[ito].pos[iwater].x;
    pos1.y = cel.atom[ito].pos[iwater].y;
    pos1.z = cel.atom[ito].pos[iwater].z;

    pos2.x = cel.atom[ito].pos[iwater2].x;
    pos2.y = cel.atom[ito].pos[iwater2].y;
    pos2.z = cel.atom[ito].pos[iwater2].z;

    pos3.x = cel.atom[ito].pos[iwater3].x;
    pos3.y = cel.atom[ito].pos[iwater3].y;
    pos3.z = cel.atom[ito].pos[iwater3].z;

    pos4.x = cel.atom[ito].pos[iwater4].x;
    pos4.y = cel.atom[ito].pos[iwater4].y;
    pos4.z = cel.atom[ito].pos[iwater4].z;

    this->put_back(pos1, pos2);
    this->put_back(pos1, pos3);
    this->put_back(pos1, pos4);

    Vector3<double> delta_pos12 = pos1-pos2;
    Vector3<double> delta_pos23 = pos2-pos3;
    Vector3<double> delta_pos43 = pos4-pos3;

    Vector3<double> delta_pos12_planar = delta_pos12 - delta_pos23*(delta_pos12*delta_pos23)/pow(delta_pos23.norm(), 2);
    Vector3<double> delta_pos43_planar = delta_pos43 - delta_pos23*(delta_pos43*delta_pos23)/pow(delta_pos23.norm(), 2);

    double angle = acos(delta_pos12_planar*delta_pos43_planar/delta_pos43_planar.norm()/delta_pos12_planar.norm())/PI*180;
    return angle;
}

void incrementalPDF2::put_back(Vector3<double> &pos1, Vector3<double> &pos2)
{
    if (abs(pos2.x-pos1.x) > INPUT.celldm1/2)
    {
        if (abs(pos2.x-pos1.x+INPUT.celldm1) <= INPUT.celldm1/2)
        {
            pos2.x += INPUT.celldm1;
        }
        if (abs(pos2.x-pos1.x-INPUT.celldm1) <= INPUT.celldm1/2)
        {
            pos2.x -= INPUT.celldm1;
        }
    }
    if (abs(pos2.y-pos1.y) > INPUT.celldm2/2)
    {
        if (abs(pos2.y-pos1.y+INPUT.celldm2) <= INPUT.celldm2/2)
        {
            pos2.y += INPUT.celldm2;
        }
        if (abs(pos2.y-pos1.y-INPUT.celldm2) <= INPUT.celldm2/2)
        {
            pos2.y -= INPUT.celldm2;
        }
    }
    if (abs(pos2.z-pos1.z) > INPUT.celldm3/2)
    {
        if (abs(pos2.z-pos1.z+INPUT.celldm3) <= INPUT.celldm3/2)
        {
            pos2.z += INPUT.celldm3;
        }
        if (abs(pos2.z-pos1.z-INPUT.celldm3) <= INPUT.celldm3/2)
        {
            pos2.z -= INPUT.celldm3;
        }
    }
    return;
}

void incrementalPDF2::output_geometry(Cell &cel, Water* water, int &iwater, int &iwater2, int &iwater3, int &iwater4, int &ito, int &ith, const int &type)
{
    Vector3<double> pos1;
    Vector3<double> pos2;
    Vector3<double> pos3;
    Vector3<double> pos4;

    Vector3<double> pos11;
    Vector3<double> pos12;
    Vector3<double> pos21;
    Vector3<double> pos22;
    Vector3<double> pos31;
    Vector3<double> pos32;
    Vector3<double> pos41;
    Vector3<double> pos42;

    pos1.x = cel.atom[ito].pos[iwater].x;
    pos1.y = cel.atom[ito].pos[iwater].y;
    pos1.z = cel.atom[ito].pos[iwater].z;

    pos2.x = cel.atom[ito].pos[iwater2].x;
    pos2.y = cel.atom[ito].pos[iwater2].y;
    pos2.z = cel.atom[ito].pos[iwater2].z;

    pos3.x = cel.atom[ito].pos[iwater3].x;
    pos3.y = cel.atom[ito].pos[iwater3].y;
    pos3.z = cel.atom[ito].pos[iwater3].z;

    pos4.x = cel.atom[ito].pos[iwater4].x;
    pos4.y = cel.atom[ito].pos[iwater4].y;
    pos4.z = cel.atom[ito].pos[iwater4].z;

    pos11.x = cel.atom[ith].pos[water[iwater].indexH[0]].x;
    pos11.y = cel.atom[ith].pos[water[iwater].indexH[0]].y;
    pos11.z = cel.atom[ith].pos[water[iwater].indexH[0]].z;

    pos12.x = cel.atom[ith].pos[water[iwater].indexH[1]].x;
    pos12.y = cel.atom[ith].pos[water[iwater].indexH[1]].y;
    pos12.z = cel.atom[ith].pos[water[iwater].indexH[1]].z;

    pos21.x = cel.atom[ith].pos[water[iwater2].indexH[0]].x;
    pos21.y = cel.atom[ith].pos[water[iwater2].indexH[0]].y;
    pos21.z = cel.atom[ith].pos[water[iwater2].indexH[0]].z;

    pos22.x = cel.atom[ith].pos[water[iwater2].indexH[1]].x;
    pos22.y = cel.atom[ith].pos[water[iwater2].indexH[1]].y;
    pos22.z = cel.atom[ith].pos[water[iwater2].indexH[1]].z;

    pos31.x = cel.atom[ith].pos[water[iwater3].indexH[0]].x;
    pos31.y = cel.atom[ith].pos[water[iwater3].indexH[0]].y;
    pos31.z = cel.atom[ith].pos[water[iwater3].indexH[0]].z;

    pos32.x = cel.atom[ith].pos[water[iwater3].indexH[1]].x;
    pos32.y = cel.atom[ith].pos[water[iwater3].indexH[1]].y;
    pos32.z = cel.atom[ith].pos[water[iwater3].indexH[1]].z;

    pos41.x = cel.atom[ith].pos[water[iwater4].indexH[0]].x;
    pos41.y = cel.atom[ith].pos[water[iwater4].indexH[0]].y;
    pos41.z = cel.atom[ith].pos[water[iwater4].indexH[0]].z;

    pos42.x = cel.atom[ith].pos[water[iwater4].indexH[1]].x;
    pos42.y = cel.atom[ith].pos[water[iwater4].indexH[1]].y;
    pos42.z = cel.atom[ith].pos[water[iwater4].indexH[1]].z;

    this->put_back(pos1, pos2);
    this->put_back(pos1, pos3);
    this->put_back(pos1, pos4);
    this->put_back(pos1, pos11);
    this->put_back(pos1, pos12);
    this->put_back(pos1, pos21);
    this->put_back(pos1, pos22);
    this->put_back(pos1, pos31);
    this->put_back(pos1, pos32);
    this->put_back(pos1, pos41);
    this->put_back(pos1, pos42);
    stringstream ss;
    ss << 12 << endl << iwater << " " << iwater2 << " " << iwater3 << " " << iwater4 << endl;
    ss << "O " << pos1.x << " " << pos1.y << " " << pos1.z << endl;
    ss << "H " << pos11.x << " " << pos11.y << " " << pos11.z << endl;
    ss << "H " << pos12.x << " " << pos12.y << " " << pos12.z << endl;
    ss << "O " << pos2.x << " " << pos2.y << " " << pos2.z << endl;
    ss << "H " << pos21.x << " " << pos21.y << " " << pos21.z << endl;
    ss << "H " << pos22.x << " " << pos22.y << " " << pos22.z << endl;
    ss << "O " << pos3.x << " " << pos3.y << " " << pos3.z << endl;
    ss << "H " << pos31.x << " " << pos31.y << " " << pos31.z << endl;
    ss << "H " << pos32.x << " " << pos32.y << " " << pos32.z << endl;
    ss << "O " << pos4.x << " " << pos4.y << " " << pos4.z << endl;
    ss << "H " << pos41.x << " " << pos41.y << " " << pos41.z << endl;
    ss << "H " << pos42.x << " " << pos42.y << " " << pos42.z << endl;
    if (type == 0 or type == 7)
    {
        this->ofs_xyz1 << ss.str().c_str();
    }
    if (type == 1 or type == 3)
    {
        this->ofs_xyz2 << ss.str().c_str();
    }
    if (type == 2 or type == 5)
    {
        this->ofs_xyz3 << ss.str().c_str();
    }
    if (type == 4 or type == 6)
    {
        this->ofs_xyz4 << ss.str().c_str();
    }
    this->ngeometry++;
    return;
}