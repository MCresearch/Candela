#include "angular_correlation.h"

angular_correlation::angular_correlation(){}

angular_correlation::~angular_correlation(){}

void angular_correlation::Routine()
{
    TITLE("angular_correlation","Routine");
    ofs_running << "Calculate the time correlation function of OH orientation." << endl;
    this->ngeo = int(round((INPUT.geo_2 - INPUT.geo_1 - INPUT.geo_ignore + 1)/INPUT.geo_interval));
    cout << "Number of geometries recorded = " << this->ngeo << endl;
    this->OH_vector = new Vector3<double>*[this->ngeo+2];
    this->nrecord = 0;
    for (int igeo=0; igeo<this->ngeo+2; igeo++)
    {
        this->OH_vector[igeo] = new Vector3<double>[INPUT.natom1];
        for (int ia=0; ia<INPUT.natom1; ia++)
        {
            this->OH_vector[igeo][ia].x = 0;
            this->OH_vector[igeo][ia].y = 0;
            this->OH_vector[igeo][ia].z = 0;
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
		++count_geometry_number;
		cout << "snapshot " << igeo << endl;
        this->calc(cel);
        this->nrecord++;
        //cout << "nrecord " << nrecord << endl;
        cel.clean();
	}//igeo
    this->output();
    for (int igeo=0; igeo<this->ngeo; igeo++)
    {
        delete[] this->OH_vector[igeo];
    }
    delete[] this->OH_vector;
    return;
}

void angular_correlation::calc(CellFile &cel)
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
        
        Vector3<double> pos_std;
        if (INPUT.func == 1)
        {
            Vector3<double> pos_H = cel.atom[ith].pos[water[iwater].indexH[0]];
            pos_std = pos_H;
            Dist2::putback_cell(cel.atom[ito].pos[iwater], pos_std);
            this->OH_vector[this->nrecord][iwater] = (pos_std - cel.atom[ito].pos[iwater])/(pos_std - cel.atom[ito].pos[iwater]).norm();
        }
        else if (INPUT.func == 2)
        {
            Vector3<double> pos_H1 = cel.atom[ith].pos[water[iwater].indexH[0]];
            Vector3<double> pos_H2 = cel.atom[ith].pos[water[iwater].indexH[1]];
            Dist2::putback_cell(cel.atom[ito].pos[iwater], pos_H1);
            Dist2::putback_cell(cel.atom[ito].pos[iwater], pos_H2);
            this->OH_vector[this->nrecord][iwater] = (pos_H1 + pos_H2 - 2.0*cel.atom[ito].pos[iwater])/(pos_H1 + pos_H2 - 2.0*cel.atom[ito].pos[iwater]).norm();
        }
        
        
    }
    delete[] water;
    return;
}

void angular_correlation::output()
{
    ofstream ofs_out ("angular_TCF.txt");
    int ndt_output = int(round(INPUT.msd_t/INPUT.msd_dt));
    double* tcf = new double[ndt_output];
    double* tcf_total = new double[ndt_output];
    for (int idt=0; idt<ndt_output; idt++) tcf_total[idt] = 0;
    for (int iwater=0; iwater<INPUT.natom1; iwater++)
    {
        for (int idt=0; idt<ndt_output; idt++) tcf[idt] = 0;
        for (int idt1=0; idt1<this->nrecord; idt1++)
        {
            for (int idt2=idt1; idt2<this->nrecord; idt2++)
            {
                int which = idt2 - idt1;
                //double delta_time = which*INPUT.msd_dt;
                if (which < ndt_output)
                {
                    tcf[which] += this->OH_vector[idt1][iwater]*this->OH_vector[idt2][iwater];
                }
            }
        }
        for (int idt=0; idt<ndt_output; idt++) tcf_total[idt] += tcf[idt];
    }
    for (int idt=0; idt<ndt_output; idt++)
    {
        ofs_out << idt*INPUT.msd_dt << " " << tcf_total[idt]/(this->nrecord-idt)/INPUT.natom1 << endl;
    }
    ofs_out.close();
    delete[] tcf;
    delete[] tcf_total;
    return;
}