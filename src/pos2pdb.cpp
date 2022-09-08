#include <iomanip>
#include "cellFile.h"
#include "input.h"
#include "pos2pdb.h"
#include "vec3.h"
#include "incremental_pdf2.h"
#include "HBs.h"
#include "gfun.h"

// added by JIANCHUAN LIU 2022-09-07

void Pos2pdb::Routine()
{
    cal();
    return;
}

void Pos2pdb::cal()
{
    cout << "convert geometry to PDB format file. Jianchuan Liu 20211013" << endl;
    const int ntype = INPUT.ntype;
    int indexatom;
    double a; // lattice parameter
    double b; // lattice parameter
    double c; // lattice parameter
    double alpha; // cosine of the angle between axis b and c
    double beta; // cosine of the angle between axis a and c
    double gamma; // cosine of the angle between axis a and b
    // create the pdb file
    ofstream ofs(INPUT.geo_out.c_str());
    // save CN  H atom number
    ofstream ofs_CN_H("CN_H_atom_num.dat");

    // cel:input geometry file
    for(int igeo=1; igeo<=INPUT.geo_2; ++igeo)
    {
        //*
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
        // read pos file 
        if( !CellFile::ReadGeometry( cel ) ) continue;
        
    	// skip the frame
        if(cel.read_and_used==false) 
		{
			cel.clean();
			continue;
		}

        if (igeo < INPUT.geo_1) // Jichuan Liu add 2021-11-04
        {
			cel.clean();
			continue;
        }
        
        //*/ 
        // if the system is water, we need to do some hydrogen bond analysis
        // so set we need to set ito, ith, and itc.
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
        }
        //
        ofs << "MODEL         " << igeo << endl;
        ofs << "REMARK   Converted from "<< INPUT.geo_directory << " file" << endl;
        ofs << "REMARK   Converted using modified D310" << endl;
        ofs << "REMARK   Snapshot time: "<< cel.snapshot_time << " ps" << endl;

        Water *water;

        // find the Cation ID
        int cation_index = -1;
        if (itca != -1 and itmg == -1)
        {
            cation_index = itca;
        }
        else if (itca == -1 and itmg != -1)
        {
            cation_index = itmg;  
        }
        assert((itca == -1 and itmg != -1) or (itca != -1 and itmg == -1));


        if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide" or INPUT.system=="2hydroxide")
        {
                
            //cout << "setup Water." << endl;
            water = new Water[cel.atom[ito].na];
            
            Water::nions = 0;
            HBs::setup_water(cel, water);
            for(int it=0; it<INPUT.ntype; ++it)
            {
                if(cel.atom[it].id != "O") continue;
                
                for(int ia=0; ia<cel.atom[it].na; ++ia)
                {
                            
                    if (water[ia].nH == 1) 
                    {
                        cout << "O indx of hydroxide: " << ia + 1 << endl;
                        ofs << "REMARK O indx of hydroxide index: "<< ia + 1  << endl;
                    }
                    if (water[ia].nH == 3) 
                    {
                        cout << "O indx of hydronium: " << ia + 1 << endl;
                        ofs << "REMARK O indx of hydronium index: "<< ia + 1  << endl;
                    }
                    for(int ih=0; ih<water[ia].nH; ++ih)
                    {
                        const int iah = water[ia].indexH[ih];
                        incrementalPDF2::put_back(cel.atom[ito].pos[ia], cel.atom[ith].pos[iah]);
                    }
                }
            }
           
        }


        ofs.width(6);
        ofs << "CRYST1";
        ofs.width(9);
        ofs.setf(ios::fixed);
        ofs.setf(ios::showpoint);
        ofs.setf(ios::right);
        a = sqrt(cel.a1.x*cel.a1.x + cel.a1.y*cel.a1.y + cel.a1.z*cel.a1.z);
        ofs << setprecision(3) << a;
        ofs.width(9);
        ofs.setf(ios::fixed);
        ofs.setf(ios::showpoint);
        ofs.setf(ios::right);
        b = sqrt(cel.a2.x*cel.a2.x + cel.a2.y*cel.a2.y + cel.a2.z*cel.a2.z);
        ofs << setprecision(3) << b;
        ofs.width(9);
        ofs.setf(ios::fixed);
        ofs.setf(ios::showpoint);
        ofs.setf(ios::right);
        c = sqrt(cel.a3.x*cel.a3.x + cel.a3.y*cel.a3.y + cel.a3.z*cel.a3.z);
        ofs << setprecision(3) << c;
        alpha = vector_cos<double> (cel.a2, cel.a3); // cosine of the angle between axis b and c
        beta = vector_cos<double> (cel.a1, cel.a3); // cosine of the angle between axis a and c
        gamma = vector_cos<double> (cel.a1, cel.a2); // cosine of the angle between axis a and b
        ofs.width(7);
        ofs.setf(ios::fixed);
        ofs.setf(ios::showpoint);
        ofs.setf(ios::right);
        ofs << setprecision(2) << alpha;
        ofs.width(7);
        ofs.setf(ios::fixed);
        ofs.setf(ios::showpoint);
        ofs.setf(ios::right);
        ofs << setprecision(2) << beta;
        ofs.width(7);
        ofs.setf(ios::fixed);
        ofs.setf(ios::showpoint);
        ofs.setf(ios::right);
        ofs << setprecision(2) << gamma << endl;

        // Select O atom around Cation ions
        int CN1[3][30];

        for  (int j = 0; j < 3; j++)
        {
            for (int i = 0; i < 30; i++)
            {
                CN1[j][i] = -1;
            }
        }
        int num_cn = 0;
        int num_cnH = 0;
        CN1[2][0] = 0;

        // output the position(x, y, z) of all atom 
        indexatom = 1;
        for(int it=0; it<ntype; ++it)
        {
        
            if (INPUT.dis_oc != -1 and cel.atom[it].id == "O")
            {
                for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
                {
                    double dist_O_C = distance(cel.atom[it].pos[ia2], cel.atom[cation_index].pos[0], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);
                    if (dist_O_C  <= INPUT.dis_oc)
                    {
                        CN1[0][num_cn] = ia2;
                        num_cn += 1;

                        if (water[ia2].nH == 1)
                        {
                            CN1[1][num_cnH] = water[ia2].indexH[0];
                            num_cnH += 1;
                            //cout << "num_cnH " << num_cnH << endl;
                        }
                        else if (water[ia2].nH == 2)
                        {
                            CN1[1][num_cnH] = water[ia2].indexH[0];
                            num_cnH += 1;
                            CN1[1][num_cnH] = water[ia2].indexH[1];
                            num_cnH += 1;
                            //cout << "num_cnH " << num_cnH << endl;
                        }
                        else
                        {
                            //assert(false);
                        }
                    }
                }
            }

            if (INPUT.dis_oc == -1)
            {
                for(int ia2=0; ia2<cel.atom[it].na; ++ia2)
                {
                    //cout << 1111 << endl;
                    // wirte the postion of pdb file
                    ofs.width(6); // 1-6
                    ofs << "ATOM  ";
                    //cout << 2222 << endl;           
                    ofs.width(5);// 7-11
                    ofs.setf(ios::right);
                    ofs << indexatom;

                    ofs.width(1); // 12
                    ofs.setf(ios::left);
                    ofs << " ";

                    ofs.width(4); // 13-16
                    ofs.setf(ios::left);
                    ofs << cel.atom[it].id;
                    
                    ofs.width(1); // 17
                    ofs.setf(ios::left);
                    ofs << " ";
                    
                    ofs.width(3); // 18-20
                    ofs.setf(ios::left);
                    ofs << "MOL";

                    ofs.width(2);// 21-22
                    ofs.setf(ios::right);
                    ofs << "A";

                    ofs.width(4);// 23-26
                    ofs.setf(ios::right);
                    ofs << "1";

                    ofs.width(4);// 26-30
                    ofs.setf(ios::right);
                    ofs << " ";

                    ofs.width(8);// 31-38
                    ofs.setf(ios::fixed);
                    ofs.setf(ios::showpoint);
                    ofs.setf(ios::right);
                    ofs << setprecision(3) << cel.atom[it].pos[ia2].x;

                    ofs.width(8);// 39-46
                    ofs.setf(ios::fixed);
                    ofs.setf(ios::showpoint);
                    ofs.setf(ios::right);
                    ofs << setprecision(3) << cel.atom[it].pos[ia2].y;
                    
                    ofs.width(8);// 47-54
                    ofs.setf(ios::fixed);
                    ofs.setf(ios::showpoint);
                    ofs.setf(ios::right);
                    ofs << setprecision(3) << cel.atom[it].pos[ia2].z;

                    ofs << endl;

                    indexatom++;
                    if(indexatom > 99999)
                    {
                        indexatom = 1;
                    }
                }
            }
        }

        // around Cation ions
        int itmp;
        if (INPUT.dis_oc != -1)
        {
            for  (int j = 0; j < 3; j++)
            {
                if (j == 0) itmp = ito;
                else if (j == 1) itmp = ith;
                else if (j == 2 and itca != -1) itmp = itca;
                else if (j == 2 and itmg != -1) itmp = itmg;
                for (int i = 0; i < num_cnH; i++)
                {
                    //cout << 1111 << endl;
                    if (CN1[j][i] == -1) continue;
                    //cout << 22222 << endl;
                    // wirte the postion of pdb file
                    ofs.width(6); // 1-6
                    //cout << 333333 << endl;
                    ofs << "ATOM  ";
                    
                    ofs.width(5);// 7-11
                    ofs.setf(ios::right);
                    
                    ofs << indexatom;

                    ofs.width(1); // 12
                    ofs.setf(ios::left);
                    ofs << " ";

                    ofs.width(4); // 13-16
                    ofs.setf(ios::left);
                    
                    ofs << cel.atom[itmp].id;
                    
                    ofs.width(1); // 17
                    ofs.setf(ios::left);
                    ofs << " ";
                    
                    ofs.width(3); // 18-20
                    ofs.setf(ios::left);
                    ofs << "MOL";

                    ofs.width(2);// 21-22
                    ofs.setf(ios::right);
                    ofs << "A";

                    ofs.width(4);// 23-26
                    ofs.setf(ios::right);
                    ofs << "1";

                    ofs.width(4);// 26-30
                    ofs.setf(ios::right);
                    ofs << " ";

                    ofs.width(8);// 31-38
                    ofs.setf(ios::fixed);
                    ofs.setf(ios::showpoint);
                    ofs.setf(ios::right);
                    ofs << setprecision(3) << cel.atom[itmp].pos[CN1[j][i]].x;

                    ofs.width(8);// 39-46
                    ofs.setf(ios::fixed);
                    ofs.setf(ios::showpoint);
                    ofs.setf(ios::right);
                    ofs << setprecision(3) << cel.atom[itmp].pos[CN1[j][i]].y;
                    
                    ofs.width(8);// 47-54
                    ofs.setf(ios::fixed);
                    ofs.setf(ios::showpoint);
                    ofs.setf(ios::right);
                    ofs << setprecision(3) << cel.atom[itmp].pos[CN1[j][i]].z;

                    ofs << endl;

                    indexatom++;
                    if(indexatom > 99999)
                    {
                        indexatom = 1;
                    }
                    
                }
            }
        }

        ofs << "TER" << endl;
        ofs << "ENDMDL" << endl;
        // output CN and number of H atom
        ofs_CN_H << cel.snapshot_time << " " << num_cn << " " << num_cnH << endl;
        cel.clean();
        if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide" or INPUT.system=="2hydroxide")
        {
            delete[] water;
        }
    }
    ofs.close();
    ofs_CN_H.close();
    return;
}

