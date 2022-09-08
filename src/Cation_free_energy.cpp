#include <iomanip>
#include "cellFile.h"
#include "input.h"
#include "Cation_free_energy.h"
#include "incremental_pdf2.h"
#include "HBs.h"
#include "gfun.h"
#include "math.h"

// added by JIANCHUAN LIU 2022-01-19

void CationFreeEnergy::Routine()
{
    cal();
    return;
}

void CationFreeEnergy::cal()
{
    // 
    cout <<"Cation Free Energy" << endl;
    const int ntype = INPUT.ntype;
    // free energy for dis_OH-OH and dis_Cation-OH
    int nxy = INPUT.rcut/INPUT.dr;
	
    this->count_xy = new double*[nxy];
	for(int ix=0; ix<nxy; ++ix)
	{
		this->count_xy[ix] = new double[nxy];
        for(int iy=0; iy<nxy; ++iy)
        {
            count_xy[ix][iy] = 0;
        }
	}

    cout << "count_xy[ix][iy]: " << count_xy[1][1] << endl;

    // for all frame
    for(int igeo=1; igeo<=INPUT.geo_2; ++igeo)
    {
        //*
		cout << "igeo=" << igeo << endl;
		CellFile cel;
        //cout << 00000000000 << endl;
		if(igeo<INPUT.geo_ignore || igeo%INPUT.geo_interval!=0) 
		{
			cel.read_and_used=false;
		}
		else cel.read_and_used=true;
		cout << "Succeeded" << endl;
		stringstream ss; ss << igeo;
		cel.file_name = ss.str();
        // read pos file 
        //cout << 111111111111 << endl;
        //CellFile::ReadGeometry( cel );
        if( !CellFile::ReadGeometry( cel ) )
        {
            //cout << 111111111111 << endl;
            continue;
            
        }
    	// skip the frame
        if(cel.read_and_used==false) 
		{
			cel.clean();
			continue;
		}

        if (igeo < INPUT.geo_1) 
        {
			cel.clean();
			continue;
        }
        
        int ito=-1;
        int ith=-1;
        int itc=-1;
        int itcl=-1;
        int itna=-1;
        int itca=-1;
        int itmg=-1;
        for(int it=0;it <INPUT.ntype; ++it)
	    {
            if(cel.atom[it].id=="O") ito=it;
            else if(cel.atom[it].id=="H" or cel.atom[it].id=="D") ith=it;
            else if(cel.atom[it].id=="C") itc=it;
            else if(cel.atom[it].id=="Cl") itcl=it;
            else if(cel.atom[it].id=="Na") itna=it;
            else if(cel.atom[it].id=="Ca") itca=it;
            else if(cel.atom[it].id=="Mg") itmg=it;

            //cout << cel.atom[it].id << endl;
	    }
        Water *water;

        if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide" or INPUT.system=="2hydroxide")
        {
            //cout << "setup Water." << endl;
            water = new Water[cel.atom[ito].na];
            Water::nions = 0;
            //cout << 111111111 << endl;
            HBs::setup_water(cel, water);
            //cout << 222222222 << endl;
        }

        // find the Cation ID
        int cation_index;
        if (itca != -1 and itmg == -1)
        {
            cation_index = itca;
        }
        else if (itca == -1 and itmg != -1)
        {
            cation_index = itmg;  
        }
        assert((itca == -1 and itmg != -1) or (itca != -1 and itmg == -1));
        //cout << "Cation index: " << getAtomIndex(cel, cation_index) << endl;
        
        int OHindex1 = -1;
        int OHindex2 = -1;
        for(int ia=0; ia<cel.atom[ito].na; ++ia) // find two OH index
        {
            if (water[ia].nH == 1 and OHindex1 == -1) OHindex1 = ia;
            else if (water[ia].nH == 1 and OHindex1 != -1 and OHindex2 == -1) OHindex2 = ia;
        }

        if (OHindex1 == -1 or OHindex2 == -1)
        {
            delete[] water;
            cel.clean();
			continue;
        }

        cout << "OHindex1: " << OHindex1 << " OHindex2: " << OHindex2 << endl;
        double dist1_O_C = distance(cel.atom[ito].pos[OHindex1], cel.atom[cation_index].pos[0], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
        double dist2_O_C = distance(cel.atom[ito].pos[OHindex2], cel.atom[cation_index].pos[0], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
        double dist_O_O = distance(cel.atom[ito].pos[OHindex1], cel.atom[ito].pos[OHindex2], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
        
        int indexX; 
        int indexY; 
        if (INPUT.cfe_model == 1) //1: X: distance of OH1-Cation, Y: distance of OH2-Cation
        {
            indexX = dist1_O_C / INPUT.dr;
            indexY = dist2_O_C / INPUT.dr;
            if(indexX<nxy and indexY<nxy and indexX>=0 and indexY>=0) count_xy[indexX][indexY] += 1;
        }
        
        if (INPUT.cfe_model == 2) //2: X: distance of OH1-Cation or OH2-Cation, Y: distance of OH-OH
        {
            indexX = dist1_O_C / INPUT.dr;
            indexY = dist_O_O / INPUT.dr;
            if(indexX<nxy and indexY<nxy and indexX>=0 and indexY>=0) count_xy[indexX][indexY] += 1;
            indexX = dist2_O_C / INPUT.dr;
            indexY = dist_O_O / INPUT.dr;
            if(indexX<nxy and indexY<nxy and indexX>=0 and indexY>=0) count_xy[indexX][indexY] += 1;
        }
        
        if (INPUT.cfe_model == 3) // 3: X: avg distance of OH1-Cation and OH2-Cation, Y:  distance of OH-OH
        {
            indexX = (dist1_O_C + dist2_O_C) / 2 / INPUT.dr;
            indexY = dist_O_O / INPUT.dr;
            if(indexX<nxy and indexY<nxy and indexX>=0 and indexY>=0) count_xy[indexX][indexY] += 1;
        }

        // clean up
        cel.clean();
        if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide" or INPUT.system=="2hydroxide")
        {
            delete[] water;
        }
    }
    //
    int sum;
	for(int ix=0; ix<nxy; ++ix)
	{
        for(int iy=0; iy<nxy; ++iy)
        {
            sum += count_xy[ix][iy];
        }
	}

    // output free energy   
    // x axis is distance between OH and Cation
    // y axis is distance between OH and OH
    ofstream ofs_result("CationFreeEnergy.dat");
	for(int iy=0; iy<nxy; ++iy)
	{
        for(int ix=0; ix<nxy; ++ix)
        {
            if(count_xy[ix][iy]>0)
            {
                count_xy[ix][iy] = -std::log(count_xy[ix][iy]/sum*INPUT.factor);
            }
            ofs_result << count_xy[ix][iy] << " ";
        }
		ofs_result << endl;
	}
	ofs_result.close();

    //  clean up
    for(int ix=0; ix<nxy; ++ix)
	{
		delete[] count_xy[ix];
	}
    delete[] count_xy;
    
}
