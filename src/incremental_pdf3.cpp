#include "cellFile.h"
#include "input.h"
#include "incremental_pdf3.h"
#include "math.h"
#include "HBs.h"
#include "gfun.h"

incrementalPDF3::incrementalPDF3(){}
incrementalPDF3::~incrementalPDF3(){}
void incrementalPDF3::Routine()
{
	TITLE("incrementalPDF3", "Routine");
	cout << "Calculate the radial distribution functions g(r)." << endl;

	assert(INPUT.dr>0.0);

	assert(INPUT.rcut>0.0);

	int count_geometry_number = 0;
    this->allocate();
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

    }//igeo 

    this->out(count_geometry_number);
    return;
}

void incrementalPDF3::allocate()
{
    this->npaths = new int[8];
    this->one_path = new int[4];
    this->two_paths = new int[10];
    this->calculated = new bool[INPUT.natom1];
    this->gr_one_path = new double*[4];
    this->gr_two_paths = new double*[10];
    this->gr_three_paths = new double[int(INPUT.rcut/INPUT.dr)];
    for (int ii=0; ii<8; ii++)
    {
        this->npaths[ii] = 0;
    }
    for (int ir=0; ir<int(INPUT.rcut/INPUT.dr); ir++)
    {
        this->gr_three_paths[ir] = 0;
    }
    for (int ii=0; ii<4; ii++)
    {
        this->one_path[ii] = 0;
        this->gr_one_path[ii] = new double[int(INPUT.rcut/INPUT.dr)];
        for (int ir=0; ir<int(INPUT.rcut/INPUT.dr); ir++)
        {
            this->gr_one_path[ii][ir] = 0;
        }
    }
    for (int ii=0; ii<10; ii++)
    {
        this->two_paths[ii] = 0;
        this->gr_two_paths[ii] = new double[int(INPUT.rcut/INPUT.dr)];
        for (int ir=0; ir<int(INPUT.rcut/INPUT.dr); ir++)
        {
            this->gr_two_paths[ii][ir] = 0;
        }
    }
    for (int ia=0; ia<INPUT.natom1; ia++)
    {
        this->calculated[ia] = false;
    }
    this->iwater_in_shell = new int[INPUT.natom1];
    this->type_of_path = new int[8];
    return;
}

void incrementalPDF3::set_false(Cell &cel, Water *water, int &iwater)
{
    for (int ia=0; ia<INPUT.natom1; ia++)
    {
        this->calculated[ia] = false;
    }
    this->calculated[iwater] = true;
    for (int idon=0; idon<water[iwater].ndonate; idon++)
    {
        int iwater2 = water[iwater].donateO[idon];
        this->calculated[iwater2] = true;
        for (int idon2=0; idon2<water[iwater2].ndonate; idon2++)
        {
            int iwater3 = water[iwater2].donateO[idon2];
            this->calculated[iwater3] = true;
        }
        for (int iacc2=0; iacc2<water[iwater2].naccept; iacc2++)
        {
            int iwater3 = water[iwater2].acceptO[iacc2];
            this->calculated[iwater3] = true;
        }
    }
    for (int iacc=0; iacc<water[iwater].naccept; iacc++)
    {
        int iwater2 = water[iwater].acceptO[iacc];
        this->calculated[iwater2] = true;
        for (int idon2=0; idon2<water[iwater2].ndonate; idon2++)
        {
            int iwater3 = water[iwater2].donateO[idon2];
            this->calculated[iwater3] = true;
        }
        for (int iacc2=0; iacc2<water[iwater2].naccept; iacc2++)
        {
            int iwater3 = water[iwater2].acceptO[iacc2];
            this->calculated[iwater3] = true;
        }
    }
    for (int ia=0; ia<INPUT.natom1; ia++)
    {
        this->iwater_in_shell[ia] = -1;
    }
    for (int ii=0; ii<8; ii++)
    {
        this->type_of_path[ii] = -1;
    }
    this->nwater_in_shell = 0;
    this->npaths_to_water = 0;
    return;
}

void incrementalPDF3::set_zero()
{
    this->npaths_to_water = 0;
    for (int ii=0; ii<8; ii++)
    {
        this->type_of_path[ii] = -1;
    }
    return;
}

void incrementalPDF3::calc(Cell &cel)
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

    for (int iwater=0; iwater<cel.atom[ito].na; iwater++)
    {
        this->set_false(cel, water, iwater);
        this->search_third_shell_water(cel, water, iwater, ito, ith);
        for (int iselect=0; iselect<this->nwater_in_shell; iselect++)
        {
            int iwater_select = this->iwater_in_shell[iselect];
            this->set_zero();
            this->search_paths(cel, water, iwater, iwater_select, ito, ith);
            if (this->npaths_to_water > 8)
            {
                cout << "There are more than 8 3-step paths connecting the two water molecules " << iwater << " " << iwater_select << endl;
                continue;
            }
            this->npaths[this->npaths_to_water-1]++;
            double dist = distance(cel.atom[ito].pos[iwater], cel.atom[ito].pos[iwater_select], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
            int which = int (dist/INPUT.dr);
            if (this->npaths_to_water == 1)
            {
                if (dist < INPUT.rcut)
                {
                    this->gr_one_path[this->type_of_path[0]][which]++;
                }
                this->one_path[this->type_of_path[0]]++;
            }
            if (this->npaths_to_water == 2)
            {
                if (this->type_of_path[0] == 0 and this->type_of_path[1] == 0)
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[0][which]++;
                    }
                    this->two_paths[0]++;
                }
                if ((this->type_of_path[0] == 0 and this->type_of_path[1] == 1) or (this->type_of_path[0] == 1 and this->type_of_path[1] == 0))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[1][which]++;
                    }
                    this->two_paths[1]++;
                }
                if ((this->type_of_path[0] == 0 and this->type_of_path[1] == 2) or (this->type_of_path[0] == 2 and this->type_of_path[1] == 0))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[2][which]++;
                    }
                    this->two_paths[2]++;
                }
                if ((this->type_of_path[0] == 0 and this->type_of_path[1] == 3) or (this->type_of_path[0] == 3 and this->type_of_path[1] == 0))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[3][which]++;
                    }
                    this->two_paths[3]++;
                }
                if ((this->type_of_path[0] == 1 and this->type_of_path[1] == 1))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[4][which]++;
                    }
                    this->two_paths[4]++;
                }
                if ((this->type_of_path[0] == 1 and this->type_of_path[1] == 2) or (this->type_of_path[0] == 2 and this->type_of_path[1] == 1))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[5][which]++;
                    }
                    this->two_paths[5]++;
                }
                if ((this->type_of_path[0] == 1 and this->type_of_path[1] == 3) or (this->type_of_path[0] == 3 and this->type_of_path[1] == 1))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[6][which]++;
                    }
                    this->two_paths[6]++;
                }
                if ((this->type_of_path[0] == 2 and this->type_of_path[1] == 2))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[7][which]++;
                    }
                    this->two_paths[7]++;
                }
                if ((this->type_of_path[0] == 2 and this->type_of_path[1] == 3) or (this->type_of_path[0] == 3 and this->type_of_path[1] == 2))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[8][which]++;
                    }
                    this->two_paths[8]++;
                }
                if ((this->type_of_path[0] == 3 and this->type_of_path[1] == 3))
                {
                    if (dist < INPUT.rcut)
                    {
                        this->gr_two_paths[9][which]++;
                    }
                    this->two_paths[9]++;
                }
            }
            if (this->npaths_to_water == 3)
            {
                if (dist < INPUT.rcut)
                {
                    this->gr_three_paths[which]++;
                }
            }
        }// for iselect
    }// for iwater
    delete[] water;
    return;
}

void incrementalPDF3::search_third_shell_water(Cell &cel, Water* water, int &iwater, int &ito, int &ith)
{
    for (int idon=0; idon<water[iwater].ndonate; idon++)
    {
        int iwater2 = water[iwater].donateO[idon];
        for (int idon2 = 0; idon2 < water[iwater2].ndonate; idon2++)
        {
            int iwater3 = water[iwater2].donateO[idon2];
            if (iwater3 != iwater)
            {
                for (int idon3 = 0; idon3 < water[iwater3].ndonate; idon3++)
                {
                    int iwater4 = water[iwater3].donateO[idon3];
                    if (iwater4 != iwater and iwater4 != iwater2 and this->calculated[iwater4]==false)
                    {
                        this->iwater_in_shell[this->nwater_in_shell] = iwater4;
                        this->nwater_in_shell++;
                    }
                }
                for (int iacc3 = 0; iacc3 < water[iwater3].naccept; iacc3++)
                {
                    int iwater4 = water[iwater3].acceptO[iacc3];
                    if (iwater4 != iwater and iwater4 != iwater2 and this->calculated[iwater4]==false)
                    {
                        this->iwater_in_shell[this->nwater_in_shell] = iwater4;
                        this->nwater_in_shell++;
                    }
                }
            }
        }
        for (int iacc2 = 0; iacc2 < water[iwater2].naccept; iacc2++)
        {
            int iwater3 = water[iwater2].acceptO[iacc2];
            if (iwater3 != iwater)
            {
                for (int idon3 = 0; idon3 < water[iwater3].ndonate; idon3++)
                {
                    int iwater4 = water[iwater3].donateO[idon3];
                    if (iwater4 != iwater and iwater4 != iwater2 and this->calculated[iwater4]==false)
                    {
                        this->iwater_in_shell[this->nwater_in_shell] = iwater4;
                        this->nwater_in_shell++;
                    }
                }
                for (int iacc3 = 0; iacc3 < water[iwater3].naccept; iacc3++)
                {
                    int iwater4 = water[iwater3].acceptO[iacc3];
                    if (iwater4 != iwater and iwater4 != iwater2 and this->calculated[iwater4]==false)
                    {
                        this->iwater_in_shell[this->nwater_in_shell] = iwater4;
                        this->nwater_in_shell++;
                    }
                }
            }
        }
    }
    for (int iacc=0; iacc<water[iwater].naccept; iacc++)
    {
        int iwater2 = water[iwater].acceptO[iacc];
        for (int idon2 = 0; idon2 < water[iwater2].ndonate; idon2++)
        {
            int iwater3 = water[iwater2].donateO[idon2];
            if (iwater3 != iwater)
            {
                for (int idon3 = 0; idon3 < water[iwater3].ndonate; idon3++)
                {
                    int iwater4 = water[iwater3].donateO[idon3];
                    if (iwater4 != iwater and iwater4 != iwater2 and this->calculated[iwater4]==false)
                    {
                        this->iwater_in_shell[this->nwater_in_shell] = iwater4;
                        this->nwater_in_shell++;
                    }
                }
                for (int iacc3 = 0; iacc3 < water[iwater3].naccept; iacc3++)
                {
                    int iwater4 = water[iwater3].acceptO[iacc3];
                    if (iwater4 != iwater and iwater4 != iwater2 and this->calculated[iwater4]==false)
                    {
                        this->iwater_in_shell[this->nwater_in_shell] = iwater4;
                        this->nwater_in_shell++;
                    }
                }
            }
        }
        for (int iacc2 = 0; iacc2 < water[iwater2].naccept; iacc2++)
        {
            int iwater3 = water[iwater2].acceptO[iacc2];
            if (iwater3 != iwater)
            {
                for (int idon3 = 0; idon3 < water[iwater3].ndonate; idon3++)
                {
                    int iwater4 = water[iwater3].donateO[idon3];
                    if (iwater4 != iwater and iwater4 != iwater2 and this->calculated[iwater4]==false)
                    {
                        this->iwater_in_shell[this->nwater_in_shell] = iwater4;
                        this->nwater_in_shell++;
                    }
                }
                for (int iacc3 = 0; iacc3 < water[iwater3].naccept; iacc3++)
                {
                    int iwater4 = water[iwater3].acceptO[iacc3];
                    if (iwater4 != iwater and iwater4 != iwater2 and this->calculated[iwater4]==false)
                    {
                        this->iwater_in_shell[this->nwater_in_shell] = iwater4;
                        this->nwater_in_shell++;
                    }
                }
            }
        }
    }
    return;
}

void incrementalPDF3::search_paths(Cell &cel, Water *water, int &iwater, int &iwater_select, int &ito, int &ith)
{
    for (int idon=0; idon<water[iwater].ndonate; idon++)
    {
        int iwater2 = water[iwater].donateO[idon];
        for (int idon2 = 0; idon2 < water[iwater2].ndonate; idon2++)
        {
            int iwater3 = water[iwater2].donateO[idon2];
            if (iwater3 != iwater)
            {
                for (int idon3 = 0; idon3 < water[iwater3].ndonate; idon3++)
                {
                    int iwater4 = water[iwater3].donateO[idon3];
                    if (iwater4 == iwater_select)
                    {
                        this->type_of_path[this->npaths_to_water] = 0;
                        this->npaths_to_water++;
                    }//DDD
                }
                for (int iacc3 = 0; iacc3 < water[iwater3].naccept; iacc3++)
                {
                    int iwater4 = water[iwater3].acceptO[iacc3];
                    if (iwater4 == iwater_select)
                    {
                        this->type_of_path[this->npaths_to_water] = 1;
                        this->npaths_to_water++;
                    }//DDA
                }
            }
        }
        for (int iacc2 = 0; iacc2 < water[iwater2].naccept; iacc2++)
        {
            int iwater3 = water[iwater2].acceptO[iacc2];
            if (iwater3 != iwater)
            {
                for (int idon3 = 0; idon3 < water[iwater3].ndonate; idon3++)
                {
                    int iwater4 = water[iwater3].donateO[idon3];
                    if (iwater4 == iwater_select)
                    {
                        this->type_of_path[this->npaths_to_water] = 2;
                        this->npaths_to_water++;
                    }//DAD
                }
                for (int iacc3 = 0; iacc3 < water[iwater3].naccept; iacc3++)
                {
                    int iwater4 = water[iwater3].acceptO[iacc3];
                    if (iwater4 == iwater_select)
                    {
                        this->type_of_path[this->npaths_to_water] = 1;
                        this->npaths_to_water++;
                    }//DAA
                }
            }
        }
    }
    for (int iacc=0; iacc<water[iwater].naccept; iacc++)
    {
        int iwater2 = water[iwater].acceptO[iacc];
        for (int idon2 = 0; idon2 < water[iwater2].ndonate; idon2++)
        {
            int iwater3 = water[iwater2].donateO[idon2];
            if (iwater3 != iwater)
            {
                for (int idon3 = 0; idon3 < water[iwater3].ndonate; idon3++)
                {
                    int iwater4 = water[iwater3].donateO[idon3];
                    if (iwater4 == iwater_select)
                    {
                        this->type_of_path[this->npaths_to_water] = 3;
                        this->npaths_to_water++;
                    }//ADD
                }
                for (int iacc3 = 0; iacc3 < water[iwater3].naccept; iacc3++)
                {
                    int iwater4 = water[iwater3].acceptO[iacc3];
                    if (iwater4 == iwater_select)
                    {
                        this->type_of_path[this->npaths_to_water] = 2;
                        this->npaths_to_water++;
                    }//ADA
                }
            }
        }
        for (int iacc2 = 0; iacc2 < water[iwater2].naccept; iacc2++)
        {
            int iwater3 = water[iwater2].acceptO[iacc2];
            if (iwater3 != iwater)
            {
                for (int idon3 = 0; idon3 < water[iwater3].ndonate; idon3++)
                {
                    int iwater4 = water[iwater3].donateO[idon3];
                    if (iwater4 == iwater_select)
                    {
                        this->type_of_path[this->npaths_to_water] = 3;
                        this->npaths_to_water++;
                    }//AAD
                }
                for (int iacc3 = 0; iacc3 < water[iwater3].naccept; iacc3++)
                {
                    int iwater4 = water[iwater3].acceptO[iacc3];
                    if (iwater4 == iwater_select)
                    {
                        this->type_of_path[this->npaths_to_water] = 0;
                        this->npaths_to_water++;
                    }//AAA
                }
            }
        }
    }
    return;
}

void incrementalPDF3::out(int &count_geometry_number)
{
    // to be output:
    // npaths, gr_one_path, gr_two_path, gr_three_path, one_path, two_paths;

    // npaths
    int sum_path = 0;
    for (int ii=0; ii<8; ii++)
    {
        sum_path += this->npaths[ii];
    }
    ofs_running << "------Statistics of number of paths connecting central water and third-shell water------" << endl;
    ofs_running << "npaths" << setw(20) << "number per molecule" << setw(20) << "percentage(%)" << endl;
    double sum = double(sum_path);
    for (int ii=0; ii<8; ii++)
    {
        double number = double(this->npaths[ii])/double(count_geometry_number)/double(INPUT.natom1);
        ofs_running << ii+1 << setw(20) << number << setw(20) << double(this->npaths[ii])/sum*100 << endl;
    }

    // gr_one_path, gr_two_path, gr_three_path
    ofstream ofs_one_path("pdf_one_path.txt");
    ofstream ofs_two_path("pdf_two_path.txt");
    ofstream ofs_three_path("pdf_three_path.txt");
    double rho = double(INPUT.natom1)/INPUT.celldm1/INPUT.celldm2/INPUT.celldm3;
    for (int ir=0; ir < int(INPUT.rcut/INPUT.dr); ir++)
    {
        ofs_one_path << double(ir)*INPUT.dr + 0.5*INPUT.dr << " ";
        double vv = 4*PI*(double(ir)+0.5)*(double(ir)+0.5)*pow(INPUT.dr, 3);
        for (int ishell=0; ishell < 4; ishell++)
        {
            ofs_one_path << this->gr_one_path[ishell][ir]/rho/double(count_geometry_number)/double(INPUT.natom1)/vv << " "; 
        }
        ofs_one_path << endl;

        ofs_two_path << double(ir)*INPUT.dr + 0.5*INPUT.dr << " ";
        for (int ishell=0; ishell < 10; ishell++)
        {
            ofs_two_path << this->gr_two_paths[ishell][ir]/rho/double(count_geometry_number)/double(INPUT.natom1)/vv << " "; 
        }
        ofs_two_path << endl;

        ofs_three_path << double(ir)*INPUT.dr + 0.5*INPUT.dr << " ";
        ofs_three_path << this->gr_three_paths[ir]/rho/double(count_geometry_number)/double(INPUT.natom1)/vv << endl; 
    }
    ofs_three_path.close();
    ofs_two_path.close();
    ofs_one_path.close();

    // one path
    ofs_running << endl << "------Statistics of types of one path------" << endl;
    ofs_running << "index" << setw(20) << "type name" << setw(20) << "percentage(%)" << endl;
    sum = 0;
    for (int ii=0; ii<4; ii++)
    {
        sum += this->one_path[ii];
    }
    for (int ii=0; ii<4; ii++)
    {
        string type_name;
        if (ii==0) type_name = "DDD/AAA";
        if (ii==1) type_name = "DDA/DAA";
        if (ii==2) type_name = "DAD/ADA";
        if (ii==3) type_name = "ADD/AAD";
        ofs_running << ii << setw(20) << type_name << setw(20) << double(this->one_path[ii])/sum*100 << endl;
    }

    // two paths
    ofs_running << endl << "------Statistics of types of two path------" << endl;
    ofs_running << "index" << setw(20) << "type name" << setw(20) << "percentage(%)" << endl;
    sum = 0;
    for (int ii=0; ii<10; ii++)
    {
        sum += this->two_paths[ii];
    }
    for (int ii=0; ii<10; ii++)
    {
        string type_name;
        if (ii==0) type_name = "00";
        if (ii==1) type_name = "01";
        if (ii==2) type_name = "02";
        if (ii==3) type_name = "03";
        if (ii==4) type_name = "11";
        if (ii==5) type_name = "12";
        if (ii==6) type_name = "13";
        if (ii==7) type_name = "22";
        if (ii==8) type_name = "23";
        if (ii==9) type_name = "33";
        ofs_running << ii << setw(20) << type_name << setw(20) << double(this->two_paths[ii])/sum*100 << endl;
    }
    return;
}
