#include "HB_stat.h"
#include "input.h"
#include "cellFile.h"

HB_stat::HB_stat(){}

HB_stat::~HB_stat(){}

void HB_stat::Routine()
{
    TITLE("HB_stat","Routine");
	
	ofs_running << "Compute HB number, length, lifespan distribution." << endl;

    assert(INPUT.system == "water" or INPUT.system == "hydronium" or INPUT.system == "hydroxide");
    int nwater = -1;
    if(INPUT.id1 == "O"){nwater = INPUT.natom1;}
    else if (INPUT.id2 == "O") {nwater = INPUT.natom2;}
    assert(nwater > 0);

    // allocate space to variables
    int nsnapshots = (int) (1+(INPUT.geo_2-INPUT.geo_1-INPUT.geo_ignore+1)/INPUT.geo_interval);
    nHBs = new int[nsnapshots];
    for (int is=0; is<nsnapshots; is++)
    {
        nHBs[is] = 0;
    }
    // func_b = 2: the HB_lifespan is the life span of an HB if it reconnects within INPUT.relax_time (default is 0.5 ps)
    assert(INPUT.stay_tmax > 0 and INPUT.stay_dt > 0);
    int n_lifespan = (int) (INPUT.stay_tmax/INPUT.stay_dt+1);
    HB_lifespan_hist = new int[n_lifespan];
    for(int ih = 0; ih < n_lifespan; ih++)
    {
        HB_lifespan_hist[ih] = 0;
    }
    if (INPUT.func_b == 2)
    {
        assert(INPUT.relax_time >= 0);
    }
    
    
    assert(INPUT.vmax > 0 and INPUT.dv > 0);
    int n_relative_vel = (int) (INPUT.vmax/INPUT.dv+1);
    relative_vel_hist = new int [n_relative_vel];
    for(int ir = 0; ir<n_relative_vel; ir++)
    {
        relative_vel_hist[ir] = 0;
    }

    int n_angle = (int) (INPUT.acut_hoo/INPUT.dtheta+1);
    angle_hist = new int[n_angle];
    for(int ia=0; ia<n_angle; ia++)
    {
        angle_hist[ia] = 0;
    }

    int n_r_oh = (int) (INPUT.rcut_oo/INPUT.dr+1);
    r_oh_hist = new int[n_r_oh];
    for(int ir=0; ir<n_r_oh; ir++)
    {
        r_oh_hist[ir] = 0;
    }

    int n_r_oo = (int) (INPUT.rcut_oo/INPUT.dr+1);
    r_oo_hist = new int[n_r_oo];
    for (int ir = 0; ir<n_r_oo; ir++)
    {
        r_oo_hist[ir] = 0;
    }

    angle_r_oo_hist = new int*[n_angle];
    for (int ia=0; ia<n_angle; ia++)
    {
        angle_r_oo_hist[ia] = new int[n_r_oo];
        for (int ir=0; ir<n_r_oo; ir++)
        {
            angle_r_oo_hist[ia][ir] = 0;
        }
    }

    angle_r_oh_hist = new int*[n_angle];
    for (int ia=0; ia<n_angle; ia++)
    {
        angle_r_oh_hist[ia] = new int[n_r_oh];
        for (int ir=0; ir<n_r_oh; ir++)
        {
            angle_r_oh_hist[ia][ir] = 0;
        }
    }

    HB_lifespan = new double* [nwater];
    incoming_time = new double* [nwater];
    // func_b = 2
    last_bonded_time = new double* [nwater];
    accepted = new bool*[nwater];
    angle = new double*[nwater];
    r_oo = new double*[nwater];
    r_oh = new double*[nwater];
    for (int iwater=0; iwater<nwater; iwater++)
    {
        HB_lifespan[iwater] = new double[nwater];
        incoming_time[iwater] = new double[nwater];
        last_bonded_time[iwater] = new double[nwater];
        accepted[iwater] = new bool[nwater];
        r_oh[iwater] = new double[nwater];
        r_oo[iwater] = new double[nwater];
        angle[iwater] = new double[nwater];
        for (int iwater2=0; iwater2<nwater; iwater2++)
        {
            HB_lifespan[iwater][iwater2] = -1;
            incoming_time[iwater][iwater2] = -1;
            last_bonded_time[iwater][iwater2] = -1;
            r_oo[iwater][iwater2] = -1;
            r_oh[iwater][iwater2] = -1;
            angle[iwater][iwater2] = -1;
            accepted[iwater][iwater2] = false;
        }
    }


    int count_geometry_number=0;
    //cout << count_geometry_number << endl;
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
		
		cout << "snapshot " << igeo << endl;

		calc(cel, count_geometry_number,
        nHBs, HB_lifespan_hist, relative_vel_hist, r_oo_hist, r_oh_hist, angle_hist,
        HB_lifespan, incoming_time, last_bonded_time, accepted, r_oo, r_oh, angle, angle_r_oo_hist, angle_r_oh_hist);
        ++count_geometry_number;
        cout << "count_geometry_number " << count_geometry_number << endl;
        cel.clean();
	}//igeo

    ofstream ofs_r_oo("r_oo.txt");
    ofstream ofs_r_oh("r_oh.txt");
    ofstream ofs_theta("theta.txt");
    ofstream ofs_vel("relative_vel.txt");
    ofstream ofs_lifespan("HB_lifespan.txt");
    ofstream ofs_nHBs("nHBs.txt");
    ofstream ofs_angle_r_oo("theta_r_oo.txt");
    ofstream ofs_angle_r_oh("theta_r_oh.txt");
    cout << 1 << endl;
    // r_oo
    double sum = 0;
    for (int ii=0; ii < n_r_oo; ii++)
    {
        sum += INPUT.dr*r_oo_hist[ii];
    }
    for (int ii=0; ii < n_r_oo; ii++)
    {
        ofs_r_oo << INPUT.dr*ii+0.5*INPUT.dr << " " << r_oo_hist[ii]/sum << endl;
    }
    ofs_r_oo.close();

    // r_oh
    sum = 0;
    for (int ii=0; ii < n_r_oh; ii++)
    {
        sum += INPUT.dr*r_oh_hist[ii];
    }
    for (int ii=0; ii < n_r_oh; ii++)
    {
        ofs_r_oh << INPUT.dr*ii+0.5*INPUT.dr << " " << r_oh_hist[ii]/sum << endl;
    }
    ofs_r_oh.close();
    // theta
    sum = 0;
    for (int ii=0; ii < n_angle; ii++)
    {
        sum += INPUT.dtheta*angle_hist[ii];
    }
    for (int ii=0; ii < n_angle; ii++)
    {
        ofs_theta << INPUT.dtheta*ii+0.5*INPUT.dtheta << " " << angle_hist[ii]/sum << endl;
    }
    ofs_theta.close();


    // relative velocity
    sum = 0;
    for (int ii=0; ii < n_relative_vel; ii++)
    {
        sum += INPUT.dv*relative_vel_hist[ii];
    }
    for (int ii=0; ii < n_relative_vel; ii++)
    {
        ofs_vel << INPUT.dv*ii+0.5*INPUT.dv << " " << relative_vel_hist[ii]/sum << endl;   
    }
    ofs_vel.close();

    // HB lifespan
    sum = 0;
    for (int ii=0; ii < n_lifespan; ii++)
    {
        sum += INPUT.stay_dt * HB_lifespan_hist[ii];
    }
    for (int ii=0; ii < n_lifespan; ii++)
    {
        ofs_lifespan << INPUT.stay_dt*ii + 0.5*INPUT.stay_dt << " " << HB_lifespan_hist[ii]/sum << endl;
    }
    ofs_lifespan.close();

    for(int isnapshot=0; isnapshot<count_geometry_number-1; isnapshot++)
    {
        ofs_nHBs << isnapshot << " " << nHBs[isnapshot] << endl;
    }
    ofs_nHBs.close();

    // angle*r_oo hist
    /*
    sum = 0;
    for (int ia=0; ia<n_angle; ia++)
    {
        for (int ir=0; ir<n_r_oo; ir++)
        {
            sum += INPUT.dtheta*INPUT.dr*angle_r_oo_hist[ia][ir];
        }
    }
    */
    for (int ir=0; ir<n_r_oo; ir++)
    {
        ofs_angle_r_oo << "  " << ir*INPUT.dr+0.5*INPUT.dr;
    }
    ofs_angle_r_oo << endl;

    assert(count_geometry_number > 0);
    for (int ia=0; ia<n_angle; ia++)
    {
        ofs_angle_r_oo << ia*INPUT.dtheta+0.5*INPUT.dtheta << " ";
        for (int ir=0; ir<n_r_oo; ir++)
        {
            double rho = INPUT.natom1/INPUT.celldm1/INPUT.celldm2/INPUT.celldm3;
            double norm_factor = rho*double(count_geometry_number)*INPUT.dtheta*INPUT.dr*2*PI*(ir+0.5)*(ir+0.5)*INPUT.dr*INPUT.dr*sin((ia+0.5)*INPUT.dtheta/180*PI);
            ofs_angle_r_oo << angle_r_oo_hist[ia][ir]/norm_factor << " ";
        }
        ofs_angle_r_oo << endl;
    }
    ofs_angle_r_oo.close();

    // angle*r_oh
    sum = 0;
    for (int ia=0; ia<n_angle; ia++)
    {
        for (int ir=0; ir<n_r_oh; ir++)
        {
            sum += INPUT.dtheta*INPUT.dr*angle_r_oh_hist[ia][ir];
        }
    }

    for (int ir=0; ir<n_r_oh; ir++)
    {
        ofs_angle_r_oh << "  " << ir*INPUT.dr+0.5*INPUT.dr;
    }
    ofs_angle_r_oh << endl;

    for (int ia=0; ia<n_angle; ia++)
    {
        ofs_angle_r_oh << ia*INPUT.dtheta+0.5*INPUT.dtheta << " ";
        for (int ir=0; ir<n_r_oh; ir++)
        {
            ofs_angle_r_oh << angle_r_oh_hist[ia][ir]/sum << " ";
        }
        ofs_angle_r_oh << endl;
    }
    ofs_angle_r_oh.close();
            
/*
    for (int iwater=0; iwater<nwater; iwater++)
    {
        delete[] HB_lifespan[iwater];
        delete[] incoming_time[iwater];
        delete[] accepted[iwater];
        delete[] r_oo[iwater];
        delete[] r_oh[iwater];
        delete[] angle[iwater];
        cout << iwater << endl;
    }
    delete[] HB_lifespan;
    delete[] incoming_time;
    delete[] accepted;
    delete[] r_oo;
    delete[] r_oh;

    delete[] nHBs;
    delete[] HB_lifespan_hist;
    delete[] relative_vel_hist;
    delete[] angle_hist;
    delete[] r_oh_hist;
    delete[] r_oo_hist;*/

    return;

}

void HB_stat::calc(const Cell &cel, const int &count_geometry_number,
 int* &nHBs, int* &HB_lifespan_hist, int* &relative_vel_hist, int* &r_oo_hist, int* &r_oh_hist, int* &angle_hist, 
 double** &HB_lifespan, double** &incoming_time, double** last_bonded_time, bool** &accepted, double** &r_oo, double** &r_oh, double** &angle,
 int** &angle_r_oo_hist, int** &angle_r_oh_hist)
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
    for(int iwater1=0; iwater1<cel.atom[ito].na; iwater1++)
    {
        for (int iwater2=0; iwater2<cel.atom[ito].na; iwater2++)
        {
            if(iwater1 == iwater2){continue;}
            bool acc = false;
            int index = -1;
            for(int iacc=0; iacc<water[iwater1].naccept; iacc++)
            {
                if (water[iwater1].acceptO[iacc] == iwater2)
                {
                    acc = true;
                    index = iacc;
                    continue;
                }
            }
            if (acc)
            {
                assert(index >= 0);
                nHBs[count_geometry_number-1] += 1;
                if(accepted[iwater1][iwater2] == false)
                {
                    accepted[iwater1][iwater2] = true;
                    if (INPUT.func_b == 1)
                    {
                        incoming_time[iwater1][iwater2] = cel.snapshot_time;
                    }
                    else if(INPUT.func_b == 2)
                    {
                        assert(cel.snapshot_time-last_bonded_time[iwater1][iwater2] > 0);
                        if (cel.snapshot_time - last_bonded_time[iwater1][iwater2] <= INPUT.relax_time and last_bonded_time[iwater1][iwater2] >= 0
                        and incoming_time[iwater1][iwater2] >= 0 and HB_lifespan[iwater1][iwater2] > 0)
                        {
                            //int lifespan_index = (int) (HB_lifespan[iwater1][iwater2] / INPUT.stay_dt);
                            //HB_lifespan_hist[lifespan_index] -= 1;
                        }
                        else if (cel.snapshot_time - last_bonded_time[iwater1][iwater2] > INPUT.relax_time and last_bonded_time[iwater1][iwater2] >= 0)
                        {
                            last_bonded_time[iwater1][iwater2] = -1;
                            incoming_time[iwater1][iwater2] = cel.snapshot_time;
                        }
                        else
                        {
                            incoming_time[iwater1][iwater2] = cel.snapshot_time;
                        }
                    }
                }
                // r_oo update
                double dis_oo = water[iwater1].accept_disO[index];
                r_oo[iwater1][iwater2] = dis_oo;
                int index_dis_oo = (int) (dis_oo/INPUT.dr);
                r_oo_hist[index_dis_oo] += 1;
                
                // r_oh update
                int index_accH = water[iwater1].acceptH[index];
                double dis_oh = distance(cel.atom[ito].pos[iwater1], cel.atom[ith].pos[index_accH], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                assert(dis_oh > 0);
                r_oh[iwater1][iwater2] = dis_oh;
                int index_dis_oh = (int) (dis_oh/INPUT.dr);
                r_oh_hist[index_dis_oh] += 1;

                // angle update
                angle[iwater1][iwater2] = water[iwater1].accept_angle[index];
                int index_angle = (int) (angle[iwater1][iwater2]/INPUT.dtheta);
                angle_hist[index_angle] += 1;

                angle_r_oo_hist[index_angle][index_dis_oo]++;
                angle_r_oh_hist[index_angle][index_dis_oh]++;

                // velocity update 
                if(INPUT.vel_file != "none")
                {
                    double dvx = cel.atom[ito].vel[iwater1].x - cel.atom[ito].vel[iwater2].x;
                    double dvy = cel.atom[ito].vel[iwater1].y - cel.atom[ito].vel[iwater2].y;
                    double dvz = cel.atom[ito].vel[iwater1].z - cel.atom[ito].vel[iwater2].z;
                    double delta_v = sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
                    delta_v *= 1e4;
                    int index_v = (int) (delta_v/INPUT.dv);
                    relative_vel_hist[index_v] += 1;
                }
            }
            else
            {
                if(accepted[iwater1][iwater2] == true)
                {
                    assert(incoming_time[iwater1][iwater2] < cel.snapshot_time);
                    HB_lifespan[iwater1][iwater2] = cel.snapshot_time-incoming_time[iwater1][iwater2];
                    assert(HB_lifespan[iwater1][iwater2] < INPUT.stay_tmax);
                    int lifespan_index = (int) (HB_lifespan[iwater1][iwater2] / INPUT.stay_dt);
                    HB_lifespan_hist[lifespan_index] += 1;

                    accepted[iwater1][iwater2] = false;
                    if (INPUT.func_b == 1)
                    {
                        incoming_time[iwater1][iwater2] = -1;
                    }
                    else if (INPUT.func_b == 2)
                    {
                        last_bonded_time[iwater1][iwater2] = cel.snapshot_time;
                    }
                }
            }
        }
    }
    delete[] water;
    return;
}
