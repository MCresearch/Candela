#include "HB_stat2.h"
#include "input.h"
#include "dist2.h"
HB_stat2::HB_stat2(){}

HB_stat2::~HB_stat2(){}

void HB_stat2::Routine()
{
    TITLE("HB_stat2","Routine");
    ofs_running << "Compute HB angle distribution." << endl;

    int nDD = 0;
    int nDA = 0;
    int nAA = 0;
    int ntheta = (int) (INPUT.theta/INPUT.dtheta)+1;
    double* DD_distr = new double[ntheta];
    double* DA_distr = new double[ntheta];
    double* AA_distr = new double[ntheta];
    this->DD_distr_acc1 = new double[ntheta];
    this->DD_distr_acc2 = new double[ntheta];
    this->DD_distr_acc3 = new double[ntheta];
    this->DA_distr_acc1 = new double[ntheta];
    this->DA_distr_acc2 = new double[ntheta];
    this->DA_distr_acc3 = new double[ntheta];
    this->AA_distr_acc1 = new double[ntheta];
    this->AA_distr_acc2 = new double[ntheta];
    this->AA_distr_acc3 = new double[ntheta];
    this->HH_distr = new double[ntheta];
    this->ofs_smallDD.open("smallDD.xyz");
    this->ofs_smallDD_acc_don.open("smallDD_acc_don.txt");
    for (int itheta = 0; itheta<ntheta; itheta++)
    {
        DD_distr[itheta] = 0;
        DA_distr[itheta] = 0;
        AA_distr[itheta] = 0;
        this->DD_distr_acc1[itheta] = 0;
        this->DD_distr_acc2[itheta] = 0;
        this->DD_distr_acc3[itheta] = 0;
        this->DA_distr_acc1[itheta] = 0;
        this->DA_distr_acc2[itheta] = 0;
        this->DA_distr_acc3[itheta] = 0;
        this->AA_distr_acc1[itheta] = 0;
        this->AA_distr_acc2[itheta] = 0;
        this->AA_distr_acc3[itheta] = 0;
        this->HH_distr[itheta] = 0;
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
        this->calc(cel, nDD, nDA, nAA, DD_distr, DA_distr, AA_distr);
        cel.clean();
		
	}//igeo

    this->ofs_smallDD_acc_don.close();
    this->ofs_smallDD.close();

    double nDD_d = (double) nDD;
    double nDA_d = (double) nDA;
    double nAA_d = (double) nAA;

    nDD_d  = nDD_d/(double)count_geometry_number;
    nDA_d  = nDA_d/(double)count_geometry_number;
    nAA_d  = nAA_d/(double)count_geometry_number;

    ofs_running << "nDD per geometry = " << nDD_d << endl;
    ofs_running << "nDA per geometry = " << nDA_d << endl;
    ofs_running << "nAA per geometry = " << nAA_d << endl;
    ofstream ofs_DD_DA_AA("DD_DA_AA_angle.txt");
    ofs_DD_DA_AA << "  DD DA AA" << endl;
    double DD_aver = 0;
    double DA_aver = 0;
    double AA_aver = 0;
    
    for (int itheta = 0; itheta<ntheta; itheta++)
    {
        ofs_DD_DA_AA << itheta*INPUT.dtheta+0.5*INPUT.dtheta << " "
        << DD_distr[itheta]/nDD << " " << DA_distr[itheta]/nDA << " " << AA_distr[itheta]/nAA << " "
        << DD_distr_acc1[itheta]/nDD << " " << DA_distr_acc1[itheta]/nDA << " " << AA_distr_acc1[itheta]/nAA << " "
        << DD_distr_acc2[itheta]/nDD << " " << DA_distr_acc2[itheta]/nDA << " " << AA_distr_acc2[itheta]/nAA << " "
        << DD_distr_acc3[itheta]/nDD << " " << DA_distr_acc3[itheta]/nDA << " " << AA_distr_acc3[itheta]/nAA << endl;
        DD_aver += (itheta*INPUT.dtheta+0.5*INPUT.dtheta)*DD_distr[itheta];
        DA_aver += (itheta*INPUT.dtheta+0.5*INPUT.dtheta)*DA_distr[itheta];
        AA_aver += (itheta*INPUT.dtheta+0.5*INPUT.dtheta)*AA_distr[itheta];
    }
    DD_aver /= nDD;
    DA_aver /= nDA;
    AA_aver /= nAA;
    ofs_running << "Average of DD O-O-O angle = " << DD_aver << endl;
    ofs_running << "Average of DA O-O-O angle = " << DA_aver << endl;
    ofs_running << "Average of AA O-O-O angle = " << AA_aver << endl;
    ofs_DD_DA_AA.close();

    double sum = 0;
    double average_HH = 0;
    ofstream ofs("HH_angle.txt");
    ofs << " angle abundance" << endl;
    for (int itheta=0; itheta<ntheta; itheta++)
    {
        sum += HH_distr[itheta]*INPUT.dtheta;
    }
    for (int itheta=0; itheta<ntheta; itheta++)
    {
        HH_distr[itheta] /= sum;
        average_HH += HH_distr[itheta]*(0.5*INPUT.dtheta + itheta*INPUT.dtheta)*INPUT.dtheta;
        ofs << 0.5*INPUT.dtheta + itheta*INPUT.dtheta << " " << HH_distr[itheta] << endl;
    }
    ofs_running << "Average of HH angle = " << average_HH << endl;
    ofs.close();
    return;
}

void HB_stat2::calc(const Cell &cel, int &nDD, int &nDA, int &nAA, double* &DD_distr, double* &DA_distr, double* &AA_distr)
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
        double angle_HH =  HBs::angle(cel, cel.atom[ith].pos[water[iwater].indexH[0]], cel.atom[ito].pos[iwater], cel.atom[ith].pos[water[iwater].indexH[1]]);
        int which_angle_HH = (int) (angle_HH/INPUT.dtheta);
        this->HH_distr[which_angle_HH]++;
        if (water[iwater].naccept != 0 and water[iwater].ndonate != 0)
        {
            for (int iacc = 0; iacc<water[iwater].naccept; iacc++)
            {
                for (int idon = 0; idon<water[iwater].ndonate; idon++)
                {
                    int iwater1 = water[iwater].acceptO[iacc];
                    int iwater2 = water[iwater].donateO[idon];
                    double angle = HBs::angle(cel, cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2]);
                    if(angle < INPUT.theta)
                    {
                        int which_angle = (int) (angle/INPUT.dtheta);
                        DA_distr[which_angle]++;
                        if(water[iwater].naccept == 1)
                        {
                            this->DA_distr_acc1[which_angle]++;
                        }
                        if(water[iwater].naccept == 2)
                        {
                            this->DA_distr_acc2[which_angle]++;
                        }
                        if(water[iwater].naccept == 3)
                        {
                            this->DA_distr_acc3[which_angle]++;
                        }
                        nDA++;
                    }
                    else
                    {
                       cout << "Warning! One angle is larger than " << INPUT.theta << " and is ignored." << endl;
                    }
                }
            }
        }// DA
        if (water[iwater].naccept >= 2)
        {
            for(int iacc1=0; iacc1<water[iwater].naccept-1; iacc1++)
            {
                for(int iacc2=iacc1+1; iacc2<water[iwater].naccept; iacc2++)
                {
                    if (iacc1 != iacc2)
                    {
                        int iwater1 = water[iwater].acceptO[iacc1];
                        int iwater2 = water[iwater].acceptO[iacc2];
                        double angle = HBs::angle(cel, cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2]);
                        if(angle < INPUT.theta)
                        {
                            int which_angle = (int) (angle/INPUT.dtheta);
                            AA_distr[which_angle]++;
                            if(water[iwater].naccept == 1)
                            {
                                this->AA_distr_acc1[which_angle]++;
                            }
                            if(water[iwater].naccept == 2)
                            {
                                this->AA_distr_acc2[which_angle]++;
                            }
                            if(water[iwater].naccept == 3)
                            {
                                this->AA_distr_acc3[which_angle]++;
                            }
                            nAA++;
                        }
                        else
                        {
                            cout << "Warning! One angle is larger than " << INPUT.theta << " and is ignored." << endl;
                        }
                    }
                }
            }
        }//AA
        // DD
        if (water[iwater].ndonate >= 2)
        {
            for(int idon1=0; idon1<water[iwater].ndonate-1; idon1++)
            {
                for(int idon2=idon1+1; idon2<water[iwater].ndonate; idon2++)
                {
                    if (idon1 != idon2)
                    {
                        int iwater1 = water[iwater].donateO[idon1];
                        int iwater2 = water[iwater].donateO[idon2];
                        double angle = HBs::angle(cel, cel.atom[ito].pos[iwater1], cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater2]);
                        if (angle < 60)
                        {
                            cout << "central water = " << iwater << endl;
                            cout << "iwater1 = " << iwater1 << ", iwater2 = " << iwater2 << endl;
                            cout << "angle = " << angle << endl;
                            if (water[iwater].donateH[idon1] == water[iwater].donateH[idon2])
                            {
                                ofs_running << "Same H for small angle." << endl;
                            }
                            if (water[iwater].donateH[idon1] != water[iwater].donateH[idon2])
                            {
                                ofs_running << "Different H for small angle." << endl;
                            }
                            ofs_smallDD_acc_don << water[iwater].ndonate << " " << water[iwater].naccept << endl;
                            int nat_record = water[iwater].ndonate+water[iwater].naccept+1;
                            ofs_smallDD << nat_record << endl;
                            ofs_smallDD << cel.snapshot_index << " " << cel.snapshot_time << endl;
                            ofs_smallDD << "O " << cel.atom[ito].pos[iwater].x << " " << cel.atom[ito].pos[iwater].y << " " << cel.atom[ito].pos[iwater].z << endl;
                            Vector3<double> newpos;
                            newpos = cel.atom[ith].pos[water[iwater].indexH[0]];
                            Dist2::putback_cell(cel.atom[ito].pos[iwater], newpos);
                            ofs_smallDD << "H " << newpos.x << " " << newpos.y << " " << newpos.z << endl;
                            newpos = cel.atom[ith].pos[water[iwater].indexH[1]];
                            Dist2::putback_cell(cel.atom[ito].pos[iwater], newpos);
                            ofs_smallDD << "H " << newpos.x << " " << newpos.y << " " << newpos.z << endl;
                            for (int iacc = 0; iacc<water[iwater].naccept; iacc++)
                            {
                                newpos = cel.atom[ito].pos[water[iwater].acceptO[iacc]];
                                Dist2::putback_cell(cel.atom[ito].pos[iwater], newpos);
                                ofs_smallDD << "O " << newpos.x << " " << newpos.y << " " << newpos.z << endl;
                                newpos = cel.atom[ith].pos[water[water[iwater].acceptO[iacc]].indexH[0]];
                                Dist2::putback_cell(cel.atom[ito].pos[iwater], newpos);
                                ofs_smallDD << "H " << newpos.x << " " << newpos.y << " " << newpos.z << endl;
                                newpos = cel.atom[ith].pos[water[water[iwater].acceptO[iacc]].indexH[1]];
                                Dist2::putback_cell(cel.atom[ito].pos[iwater], newpos);
                                ofs_smallDD << "H " << newpos.x << " " << newpos.y << " " << newpos.z << endl;
                            }
                            for (int idon = 0; idon<water[iwater].ndonate; idon++)
                            {
                                newpos = cel.atom[ito].pos[water[iwater].donateO[idon]];
                                Dist2::putback_cell(cel.atom[ito].pos[iwater], newpos);
                                ofs_smallDD << "O " << newpos.x << " " << newpos.y << " " << newpos.z << endl;
                                newpos = cel.atom[ith].pos[water[water[iwater].donateO[idon]].indexH[0]];
                                Dist2::putback_cell(cel.atom[ito].pos[iwater], newpos);
                                ofs_smallDD << "H " << newpos.x << " " << newpos.y << " " << newpos.z << endl;
                                newpos = cel.atom[ith].pos[water[water[iwater].donateO[idon]].indexH[1]];
                                Dist2::putback_cell(cel.atom[ito].pos[iwater], newpos);
                                ofs_smallDD << "H " << newpos.x << " " << newpos.y << " " << newpos.z << endl;
                            }
                            //exit(0);
                        }
                        if (angle < INPUT.theta)
                        {
                            int which_angle = (int) (angle/INPUT.dtheta);
                            DD_distr[which_angle]++;
                            if(water[iwater].naccept == 1)
                            {
                                this->DD_distr_acc1[which_angle]++;
                            }
                            if(water[iwater].naccept == 2)
                            {
                                this->DD_distr_acc2[which_angle]++;
                            }
                            if(water[iwater].naccept == 3)
                            {
                                this->DD_distr_acc3[which_angle]++;
                            }
                            nDD++;
                        }
                        else
                        {
                            cout << "Warning! One angle is larger than " << INPUT.theta << " and is ignored." << endl;
                        }
                    }
                }
            }
        }//DD
    }//iwater
    delete[] water;
    return;
}