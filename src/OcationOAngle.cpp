#include <iomanip>
#include "cellFile.h"
#include "input.h"
#include "OcationOAngle.h"
//#include "vec3.h"
#include "incremental_pdf2.h"
#include "HBs.h"
#include "gfun.h"

// added by JIANCHUAN LIU 2022-09-07

void OcationOAngle::Routine()
{
    cal();
    return;
}


double OcationOAngle::angle1(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &pos3) 
{
	double x1 = pos1.x; double y1 = pos1.y; double z1 = pos1.z; 
	double x2 = pos2.x; double y2 = pos2.y; double z2 = pos2.z; 
	double x3 = pos3.x; double y3 = pos3.y; double z3 = pos3.z;

	double d12x = shortest(x1, x2, INPUT.celldm1);
	double d12y = shortest(y1, y2, INPUT.celldm2);
	double d12z = shortest(z1, z2, INPUT.celldm3);

	double d32x = shortest(x3, x2, INPUT.celldm1);
	double d32y = shortest(y3, y2, INPUT.celldm2);
	double d32z = shortest(z3, z2, INPUT.celldm3);
	
	double norm1 = sqrt(d12x*d12x+d12y*d12y+d12z*d12z);
	double norm2 = sqrt(d32x*d32x+d32y*d32y+d32z*d32z);

	double angle = (d12x*d32x+d12y*d32y+d12z*d32z)/norm1/norm2;	

	angle = acos(angle)/3.1415926535897*180;

	return angle;
}


void OcationOAngle::cal()
{
    // 
    cout <<"O-Caton-O angle" << endl;
    const int ntype = INPUT.ntype;

    ofstream ofs_angle_withOH("Angle_O-Cation-O_with_1OH.dat"); // Only calculate the OH around Cation within "dis_oc"
    ofstream ofs_angle_withOH1("Angle_O-Cation-O_with_2OH.dat"); // Only calculate the OH around Cation within "dis_oc"
    ofstream ofs_angle_withoutOH("Angle_O-Cation-O_without_OH.dat"); // Only calculate the OH around Cation within "dis_oc"
    ofstream ofs_cn_withOH("Coordination_numbers_with_1OH.dat"); // Only calculate the OH around Cation within "dis_oc"
    ofstream ofs_cn_withOH1("Coordination_numbers_with_2OH.dat"); // Only calculate the OH around Cation within "dis_oc"
    ofstream ofs_cn_withoutOH("Coordination_numbers_without_OH.dat"); // Only calculate the OH around Cation within "dis_oc"
    ofstream ofs_vector_withOH("Vector_with_1OH.dat"); // Only calculate the OH around Cation within "dis_oc"
    ofstream ofs_vector_withOH1("Vector_with_2OH.dat"); // Only calculate the OH around Cation within "dis_oc"
    ofstream ofs_vector_withoutOH("Vector_without_OH.dat"); // Only calculate the OH around Cation within "dis_oc"

    double O_Cation_O_with[90];
    double Vector_with[90];
    double O_Cation_O_with1[90];
    double Vector_with1[90];
    double O_Cation_O_without[90];
    double Vector_without[90];
    for (int i = 0; i<90; i++)
    {
        O_Cation_O_with[i] = 0;
        Vector_with[i] = 0;
        O_Cation_O_with1[i] = 0;
        Vector_with1[i] = 0;
        O_Cation_O_without[i] = 0;
        Vector_without[i] = 0;
    }

    double Coordination_numbers_with[10];
    double Coordination_numbers_with1[10];
    double Coordination_numbers_without[10];
    for (int i = 0; i<10; i++)
    {
        Coordination_numbers_with[i] = 0;
        Coordination_numbers_with1[i] = 0;
        Coordination_numbers_without[i] = 0;
    }



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
            HBs::setup_water(cel, water);
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
        cout << "Cation index: " << getAtomIndex(cel, cation_index) << endl;

        // Select O atom around Cation ions
        int CN[10];
        for (int i = 0; i < 10; i++)
        {
            CN[i] = -1;
        }
        int num_cn = 0;

        double dist_O_C;
        bool with_OH1 = false;
        bool with_OH2 = false;
        for(int it=0; it<INPUT.ntype; ++it)
        {
            if(cel.atom[it].id == "O")
            {
                for(int ia=0; ia<cel.atom[it].na; ++ia)
                {
                    dist_O_C = distance(cel.atom[it].pos[ia], cel.atom[cation_index].pos[0], INPUT.celldm1, INPUT.celldm2, INPUT.celldm3);

                    // Coordination numbers ID 
                    if (dist_O_C  <= INPUT.dis_oc)
                    {
                        // find OH1
                        if (water[ia].nH == 1 and with_OH1 == false and with_OH2 == false)
                        {
                            with_OH1 = true;
                        }
                        // find OH2
                        else if (water[ia].nH == 1 and with_OH1 == true and with_OH2 == false)
                        {
                            with_OH2 = true;
                        }
                        CN[num_cn] = ia;
                        num_cn += 1;
                    }

                } // ia
            } // select O atom
        } // it
        
        // Calculate the Angle of O-Cation-O
        cout << "Coordination numbers of Cation: " << num_cn << endl;
        for (int i = 0; i < num_cn; i++)
        {
            for (int j = i+1; j < num_cn; j++)
            {
                //cout << i << " " << j << endl;
                // O-Cation-O angle
                double angle_O_Cation_O = angle1(cel, cel.atom[ito].pos[CN[i]], cel.atom[cation_index].pos[0], cel.atom[ito].pos[CN[j]]);
                //cout << "O-Cation-O angle: " << angle_O_Cation_O << endl;
                //with one OH
                if (with_OH1 == true and with_OH2 == false)
                {
                    //cout << "O-Cation-O angle: " << angle_O_Cation_O << endl;
                    ofs_angle_withOH << angle_O_Cation_O << endl;
                    for (int nn = 0; nn<90; nn++)
                    {
                        if((angle_O_Cation_O > nn*2) and (angle_O_Cation_O <= (nn + 1)*2)) O_Cation_O_with[nn] += 1;
                    }
                }
                // with two OH
                if (with_OH1 == true and with_OH2 == true)
                {
                    //cout << "O-Cation-O angle: " << angle_O_Cation_O << endl;
                    ofs_angle_withOH1 << angle_O_Cation_O << endl;
                    for (int nn = 0; nn<90; nn++)
                    {
                        if((angle_O_Cation_O > nn*2) and (angle_O_Cation_O <= (nn + 1)*2)) O_Cation_O_with1[nn] += 1;
                    }
                }
                //without OH (pure water)
                else if (with_OH1 == false and with_OH2 == false )
                {
                    //cout << "O-Cation-O angle: " << angle_O_Cation_O << endl;
                    ofs_angle_withoutOH << angle_O_Cation_O << endl;
                    for (int nn = 0; nn<90; nn++)
                    {
                        if((angle_O_Cation_O > nn*2) and (angle_O_Cation_O <= (nn + 1)*2)) O_Cation_O_without[nn] += 1;
                    }
                }

            }
        }
        // calculate the angle of Cation-O vector
        for (int ia = 0; ia < num_cn; ia++)
        {
            if (water[CN[ia]].nH == 2)
            {
                //with one OH
                if (with_OH1 == true and with_OH2 == false)
                {
                    int indexH1 = water[CN[ia]].indexH[0];
                    int indexH2 = water[CN[ia]].indexH[1];
                    Vector3<double> vector1;
                    angle_vector(cel, cel.atom[ith].pos[indexH1],cel.atom[ito].pos[CN[ia]], cel.atom[ith].pos[indexH2], vector1);
                    //cout << "v11:" << vector1[0] << "v12:" << vector1[1] << "v13:" << vector1[3] << endl;
                    Vector3<double> vector2;
                    atom_vector(cel, cel.atom[ito].pos[CN[ia]], cel.atom[cation_index].pos[0], vector2);
                    double angle = angle3(vector1.x, vector1.y, vector1.z, vector2.x, vector2.y, vector2.z);
                    cout << "Angle with Cation: " << angle << endl;

                    ofs_vector_withOH << angle << endl;
                    for (int nn = 0; nn<90; nn++)
                    {
                        if((angle > nn*2) and (angle <= (nn + 1)*2)) Vector_with[nn] += 1;
                    }
                }
                //with two OH
                else if (with_OH1 == true and with_OH2 == true)
                {
                    int indexH1 = water[CN[ia]].indexH[0];
                    int indexH2 = water[CN[ia]].indexH[1];
                    Vector3<double> vector1;
                    angle_vector(cel, cel.atom[ith].pos[indexH1],cel.atom[ito].pos[CN[ia]], cel.atom[ith].pos[indexH2], vector1);
                    //cout << "v11:" << vector1[0] << "v12:" << vector1[1] << "v13:" << vector1[3] << endl;
                    Vector3<double> vector2;
                    atom_vector(cel, cel.atom[ito].pos[CN[ia]], cel.atom[cation_index].pos[0], vector2);
                    double angle = angle3(vector1.x, vector1.y, vector1.z, vector2.x, vector2.y, vector2.z);
                    cout << "Angle with Cation: " << angle << endl;

                    ofs_vector_withOH1 << angle << endl;
                    for (int nn = 0; nn<90; nn++)
                    {
                        if((angle > nn*2) and (angle <= (nn + 1)*2)) Vector_with1[nn] += 1;
                    }
                }
                //without OH (pure water)
                else if (with_OH1 == false and with_OH2 == false)
                {
                    int indexH1 = water[CN[ia]].indexH[0];
                    int indexH2 = water[CN[ia]].indexH[1];
                    //cout << "!!!! test: " << indexH1 << endl;
                    //cout << "!!!! test: " << indexH2 << endl;

                    Vector3<double> vector1;
                    angle_vector(cel, cel.atom[ith].pos[indexH1], cel.atom[ito].pos[CN[ia]], cel.atom[ith].pos[indexH2], vector1);
                    //cout << "v11:" << vector1[0] << "v12:" << vector1[1] << "v13:" << vector1[3] << endl;
                    Vector3<double> vector2;
                    atom_vector(cel, cel.atom[ito].pos[CN[ia]], cel.atom[cation_index].pos[0], vector2);
                    double angle = angle3(vector1.x, vector1.y, vector1.z, vector2.x, vector2.y, vector2.z);
                    //cout << "vector1.x: " << vector1.x << " vector1.y: " << vector1.y << " vector1.z: " << vector1.z  << endl;
                    //cout << "vector2.x: " << vector2.x << " vector2.y: " << vector2.y << " vector2.z: " << vector2.z  << endl;
                    cout << "Angle without Cation: " << angle << endl;
                    ofs_vector_withoutOH << angle << endl;
                    for (int nn = 0; nn<90; nn++)
                    {
                        if((angle > nn*2) and (angle <= (nn + 1)*2)) Vector_without[nn] += 1;
                    }
                }
            }
        }
        // Coordination numbers with one OH
        if (with_OH1 == true and with_OH2 == false)
        {
            //cout << "O-Cation-O angle: " << angle_O_Cation_O << endl;
            ofs_cn_withOH << num_cn << endl;
            Coordination_numbers_with[num_cn] += 1;
        }
        // Coordination numbers with two OH
        if (with_OH1 == true and with_OH2 == true)
        {
            //cout << "O-Cation-O angle: " << angle_O_Cation_O << endl;
            ofs_cn_withOH1 << num_cn << endl;
            Coordination_numbers_with1[num_cn] += 1;
        }
        // Coordination numbers without OH (pure water)
        else if (with_OH1 == false and with_OH2 == false)
        {
            //cout << "O-Cation-O angle: " << angle_O_Cation_O << endl;
            ofs_cn_withoutOH << num_cn << endl;
            Coordination_numbers_without[num_cn] += 1;
        }
        //
        cel.clean();
        if(INPUT.system=="water" or INPUT.system=="hydronium" or INPUT.system=="hydroxide" or INPUT.system=="2hydroxide")
        {
            delete[] water;
        }
    }

    //

    int total_O_Cation_O_with = 0;
    int total_Vector_with = 0;
    int total_O_Cation_O_with1 = 0;
    int total_Vector_with1 = 0;
    int total_O_Cation_O_without = 0;
    int total_Vector_without = 0;

    int total_CN_with = 0;
    int total_CN_with1 = 0;
    int total_CN_without = 0;

    for (int nn = 0; nn<90; nn++)
    {
        total_O_Cation_O_with += O_Cation_O_with[nn];
        total_O_Cation_O_with1 += O_Cation_O_with1[nn];
        total_O_Cation_O_without += O_Cation_O_without[nn];
        total_Vector_with += Vector_with[nn];
        total_Vector_with1 += Vector_with1[nn];
        total_Vector_without += Vector_without[nn];
        if (nn < 10)
        {
            total_CN_with += Coordination_numbers_with[nn];
            total_CN_with1 += Coordination_numbers_with1[nn];
            total_CN_without += Coordination_numbers_without[nn];
        }
    }

    ofstream dis_angle("dis_O-Cation-O.dat");
    ofstream dis_Vector("dis_Vector.dat");
    ofstream dis_CN("dis_Coordination_numbers.dat");

    dis_angle << "@ Degree   Dis(with1OH)  Dis(with2OH)  Dis(withoutOH)" << endl;
    dis_Vector << "@ Degree   Dis(with1OH) Dis(with2OH)  Dis(withoutOH)" << endl;
    dis_CN << "@ CN   Dis(with1OH)  Dis(with2OH)  Dis(withoutOH)" << endl;

    for (int nn = 0; nn<90; nn++)
    {
        dis_angle.width(15);
        dis_angle << (nn + 0.5) * 2;
        dis_angle.width(15);
        dis_angle << O_Cation_O_with[nn] / total_O_Cation_O_with;
        dis_angle.width(15);
        dis_angle << O_Cation_O_with1[nn] / total_O_Cation_O_with1;
        dis_angle.width(15);
        dis_angle << O_Cation_O_without[nn] / total_O_Cation_O_without << endl;

        dis_Vector.width(15);
        dis_Vector << (nn + 0.5) * 2;
        dis_Vector.width(15);
        dis_Vector << Vector_with[nn] / total_Vector_with;
        dis_Vector.width(15);
        dis_Vector << Vector_with1[nn] / total_Vector_with1;
        dis_Vector.width(15);
        dis_Vector << Vector_without[nn] / total_Vector_without << endl;

        if (nn < 10)
        {
            dis_CN.width(15);
            dis_CN << nn;
            dis_CN.width(15);
            dis_CN << Coordination_numbers_with[nn] / total_CN_with;
            dis_CN.width(15);
            dis_CN << Coordination_numbers_with1[nn] / total_CN_with1;
            dis_CN.width(15);
            dis_CN << Coordination_numbers_without[nn] / total_CN_without << endl;
        }
    } 

    dis_angle.close();
    dis_Vector.close();
    dis_CN.close();


    ofs_angle_withOH.close();
    ofs_angle_withOH1.close();
    ofs_angle_withoutOH.close();
    ofs_cn_withOH.close();
    ofs_cn_withOH1.close();
    ofs_cn_withoutOH.close();
    ofs_vector_withOH.close();
    ofs_vector_withOH1.close();
    ofs_vector_withoutOH.close();
}

// get the first atom index for every atom type, beging 0.
int OcationOAngle::getAtomIndex(const Cell &cel, const int &ia1)
{
    int index = 0;
    if (ia1 == 0)
    {
        index = 0;
    }
    else 
    {
        for (int in = 0; in <= ia1; in++)
        {
            index += cel.atom[in].na;

        }
    }
    //cout << "index: " << index << endl;
    return index;
}

// 计算水分子的角平分线向量
double  OcationOAngle::angle_vector(const Cell &cel, Vector3<double> &posH1, Vector3<double> &posO, Vector3<double> &posH2, Vector3<double> &vector1) 
{
	double x1 = posH1.x; double y1 = posH1.y; double z1 = posH1.z; 
	double x2 = posO.x; double y2 = posO.y; double z2 = posO.z; 
	double x3 = posH2.x; double y3 = posH2.y; double z3 = posH2.z;

	double d12x = shortest(x1, x2, INPUT.celldm1);
	double d12y = shortest(y1, y2, INPUT.celldm2);
	double d12z = shortest(z1, z2, INPUT.celldm3);

	double d32x = shortest(x3, x2, INPUT.celldm1);
	double d32y = shortest(y3, y2, INPUT.celldm2);
	double d32z = shortest(z3, z2, INPUT.celldm3);
	
	double norm1 = sqrt(d12x*d12x+d12y*d12y+d12z*d12z);
	double norm2 = sqrt(d32x*d32x+d32y*d32y+d32z*d32z);

    double vector[3];
    
    vector1.x = (x1-x2)/norm1 + (x3-x2)/norm2;
    vector1.y = (y1-y2)/norm1 + (y3-y2)/norm2;
    vector1.z = (z1-z2)/norm1 + (z3-z2)/norm2;
}

// 计算两个原子向量
double OcationOAngle::atom_vector(const Cell &cel, Vector3<double> &pos1, Vector3<double> &pos2, Vector3<double> &vector2)
{
    double x1 = pos1.x; double y1 = pos1.y; double z1 = pos1.z; 
	double x2 = pos2.x; double y2 = pos2.y; double z2 = pos2.z; 
    //cout << "x1: " << x1 << " y1: " << y1 << " z1: " << z1 << endl;
    //cout << "x2: " << x2 << " y2: " << y2 << " z2: " << z2 << endl;
    //double vector[3];
    vector2.x = (x1-x2);
    //cout << "vector[0]: " << vector2.x << " x1-x2: " << x1 - x2 << endl;
    vector2.y = (y1-y2);
    vector2.z = (z1-z2);
} 
// 计算两个向量夹角
double OcationOAngle::angle3(double &v1x, double &v1y, double &v1z, double &v2x, double &v2y, double &v2z)
{
	double norm1 = sqrt(v1x*v1x+v1y*v1y+v1z*v1z);
	double norm2 = sqrt(v2x*v2x+v2y*v2y+v2z*v2z);

    double ab = v1x*v2x+v1y*v2y+v1z*v2z;

    double angle = ab/norm1/norm2;	

	angle = acos(angle)/3.1415926535897*180;

    return angle;
}