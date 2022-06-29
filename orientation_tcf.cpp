#include "orientation_tcf.h"
#include "gfun.h"
Orientation_TCF::Orientation_TCF(){}

Orientation_TCF::~Orientation_TCF(){}

void Orientation_TCF::Routine()
{
    TITLE("HB_correlation","Routine");
    ofs_running << "Calculate the orientation time correlation function of water or hydroxide or hydronium." << endl;

    this->ndt = (INPUT.geo_2 - INPUT.geo_1 + 1 - INPUT.geo_ignore)/INPUT.geo_interval;
    this->ndt1 = ndt + 10;

    this->time_serial = new double[this->ndt1];

    if (INPUT.system == "water") 
    {
        this->n_mlc = INPUT.natom1;
    }
    else if (INPUT.system == "hydroxide" or INPUT.system == "hydronium")
    {
        this->n_mlc = 1;
    }

    this->orient_vec = new Vector3<double>*[this->ndt1];
    for (int idt=0; idt<this->ndt1; idt++)
    {
        this->orient_vec[idt] = new Vector3<double>[this->n_mlc];
        for (int imlc=0; imlc<this->n_mlc; imlc++)
        {
            this->orient_vec[idt][imlc].x = 0;
            this->orient_vec[idt][imlc].y = 0;
            this->orient_vec[idt][imlc].z = 0;
        }
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
		
		cout << "snapshot " << igeo << endl;
        assert(count_geometry_number < this->ndt1);
        this->calc(cel, count_geometry_number);
        ++count_geometry_number;
        
        //this->ofs_nHB << cel.snapshot_time << " " << this->nHB << endl;
        cel.clean();
	}//igeo

    this->time_serial[0] = 0;
    this->ndt = count_geometry_number;
    cout << "this->ndt = " << this->ndt;
    /*
    for (int idt=0; idt<this->ndt; idt++)
    {
        cout << "this->time_serial[" << idt <<  "] = " << this->time_serial[idt] << endl;
    }
    */
    this->calc_tcf();
    return;
}

void Orientation_TCF::calc(CellFile &cel, int &count_geometry_number)
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

    int ia_record = 0;
    for (int ia=0; ia<cel.atom[ito].na; ia++)
    {
        if (INPUT.system == "hydronium" and water[ia].nH != 3) continue;
        if (INPUT.system == "hydroxide" and water[ia].nH != 1) continue;
    
        Vector3<double> water_vec;
        water_vec.x = 0;
        water_vec.y = 0;
        water_vec.z = 0;
        for (int ih=0; ih<water[ia].nH; ih++)
        {
            Vector3<double> oh_vec;
            oh_vec.x = shortest(cel.atom[ith].pos[water[ia].indexH[ih]].x, cel.atom[ito].pos[ia].x, INPUT.celldm1);
            oh_vec.y = shortest(cel.atom[ith].pos[water[ia].indexH[ih]].y, cel.atom[ito].pos[ia].y, INPUT.celldm2);
            oh_vec.z = shortest(cel.atom[ith].pos[water[ia].indexH[ih]].z, cel.atom[ito].pos[ia].z, INPUT.celldm3);
            oh_vec = oh_vec/oh_vec.norm();

            water_vec += oh_vec;
        }
        water_vec /= water_vec.norm();
        this->orient_vec[count_geometry_number][ia_record].x = water_vec.x;
        this->orient_vec[count_geometry_number][ia_record].y = water_vec.y;
        this->orient_vec[count_geometry_number][ia_record].z = water_vec.z;
        ia_record++;
    }// for ia

    if (count_geometry_number == 0) this->time_serial[count_geometry_number] = cel.snapshot_time;
    else
    {
        this->time_serial[count_geometry_number] = cel.snapshot_time - this->time_serial[0];
    }
    delete[] water;
    return;
}

void Orientation_TCF::calc_tcf()
{
    int record_ndt = int(INPUT.msd_t/INPUT.msd_dt)+10;
    cout << "record_ndt = " << record_ndt << endl;
    int* record_num = new int [record_ndt];
    for (int it=0; it<record_ndt; it++) record_num[it] = 0;
    double* TCF = new double[record_ndt];
    for (int it=0; it<record_ndt; it++) TCF[it] = 0; 
    cout << "Allcation done." << endl;
    cout << "nmolecules = " << this->n_mlc << endl;
    cout << "this->ndt = " << this->ndt << endl;
    for (int it1=0; it1<this->ndt; it1++)
    {
        cout << "it1 = " << it1 << ", time_serial[it1] = " << this->time_serial[it1] << endl;
        for (int it2=it1; it2<this->ndt; it2++)
        {
            if (this->time_serial[it2] - this->time_serial[it1] < INPUT.msd_t and this->time_serial[it2] - this->time_serial[it1] >= 0)
            {
                //cout << "time interval = " << this->time_serial[it2] << "-" << time_serial[it1] << " " << this->time_serial[it2] - this->time_serial[it1] << endl;
                for (int imlc=0; imlc<this->n_mlc; imlc++)
                {
                    if (this->orient_vec[it1][imlc].norm() > 1e-15 and this->orient_vec[it2][imlc].norm() > 1e-15)
                    {
                        double prdct = this->orient_vec[it1][imlc]*this->orient_vec[it2][imlc];
                        int it_idx = (this->time_serial[it2] - this->time_serial[it1])/INPUT.msd_dt;
                        assert(it_idx < record_ndt and it_idx >= 0);
                        TCF[it_idx] += (3*pow(prdct, 2) - 1)/2;
                        ++record_num[it_idx];
                    }
                }
            }
        }
    }
    cout << "Calculation finished. " << endl;
    ofstream ofs("orientation_TCF.txt");

    for (int it=0; it<record_ndt; it++)
    {
        if (record_num[it] > 0)
        {
            TCF[it] /= record_num[it];
            ofs << it*INPUT.msd_dt << " " << TCF[it] << endl;
            cout << it*INPUT.msd_dt << " " << TCF[it] << endl;
        }
        else
        {
            ofs << it*INPUT.msd_dt << " " << 0 << endl;
            cout << it*INPUT.msd_dt << " " << 0 << endl;
        }
    }
    ofs.close();
    return;
}
