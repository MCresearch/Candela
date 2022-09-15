#include "gfun.h"
#include <cstdlib>

ofstream ofs_running;
int NPROC=1;
int RANK=0;

// CONSTANTS
double PI=3.1415926535897932384626;
double KB=1.38064852*1.0e-23; // Botzmann constant. unit: m^2 * kg * s^-2 * K^-1
double EPSILON0=8.854187817*1.0e-12; // Vacuum permittivity. F * m^-1, 1 F= 1C/1V
double BOHR=0.52917721092; // 1 Bohr = 0.529 Angstrom


bool SCAN_BEGIN(ifstream &ifs, const string &TargetName, const bool restart)
{
    string SearchName;
    bool find = false;
    if (restart)
    {
        ifs.clear();
        ifs.seekg(0);
    }
    ifs.rdstate();
    while (ifs.good())
    {
        ifs >> SearchName;
		//cout << " " << SearchName << endl;
        if ( SearchName == TargetName)
        {
            find = true;
            break;
        }
    }
    if (!find)
    {
     //   cout <<" In SCAN_BEGIN, can't find: "<<TargetName<<" block."<<endl;
    }
    return find;
}



void TITLE(const string &class_name,const string &function_name)
{
	//return;
    //cout<<" ==> "<<class_name<<"::"<<function_name<<endl;
    return;
}

void QUIT(const string &reason)
{
	cout << " ------------------------------------------ " << endl;
	cout << " Quit because : " << reason << endl;
	exit(0);
}


double shortest(const double &tmp_pos, const double &pos, const double &celldm)
{
	double vec1 = tmp_pos-pos;
	double vec2 = tmp_pos-pos+celldm;
	double vec3 = tmp_pos-pos-celldm;

	double dis1 = abs(vec1);
	double dis2 = abs(vec2);
	double dis3 = abs(vec3);

	if(dis1<=dis2 and dis1<=dis3) return vec1; 
	if(dis2<=dis1 and dis2<=dis3) return vec2; 
	if(dis3<=dis1 and dis3<=dis2) return vec3; 
} 

double Polynomial_Interpolation
(
    const double *table,
    const int &table_length,
    const double &table_interval,
    const double &x				// input value
)
{
//	assert(table_interval>0);
    const double position = x / table_interval;
    const int iq = static_cast<int>(position);
//	if(iq >= table_length-4)
//		cout << "\n iq = " << iq << " table_length = " << table_length;
  
   assert(iq < table_length-4);
    const double x0 = position - static_cast<double>(iq);
    const double x1 = 1.0 - x0;
    const double x2 = 2.0 - x0;
    const double x3 = 3.0 - x0;

    /*
    const double y=
    	table[iq]   * x1 * x2 * x3 / 6.0 +
    	table[iq+1] * x0 * x2 * x3 / 2.0 -
    	table[iq+2] * x1 * x0 * x3 / 2.0 +
    	table[iq+3] * x1 * x2 * x0 / 6.0 ;
    	*/

    return x1*x2*(table[iq]*x3+table[iq+3]*x0)/6.0
         + x0*x3*(table[iq+1]*x2-table[iq+2]*x1)/2.0;
}

double Polynomial_Interpolation_xy
(
    const double *xpoint,
    const double *ypoint,
    const int table_length,
    const double &x             // input value
)
{
    int position = -1;

    if (x < xpoint[0])
    {
        return ypoint[0];
    }
    // timer::tick("Mathzone","Poly_Inter_xy");

    for (int ik = 0; ik < table_length; ik++)
    {
        if (x < xpoint[ik])
        {
            break;
        }
        else
        {
            position ++;
        }
    }

    assert(position >= 0);
    assert(position <= table_length-1);

    if (position + 6 < table_length)
    {
        double dx1, dx2, dx3, dx4, dx5, dx6;
        dx1 = x - xpoint[position];
        dx2 = x - xpoint[position+1];
        dx3 = x - xpoint[position+2];
        dx4 = x - xpoint[position+3];
        dx5 = x - xpoint[position+4];
        dx6 = x - xpoint[position+5];


        double x12, x13, x14, x15, x16, x23, x24, x25, x26, x34, x35, x36, x45, x46, x56;
        x12 = xpoint[position] - xpoint[position+1];
        x13 = xpoint[position] - xpoint[position+2];
        x14 = xpoint[position] - xpoint[position+3];
        x15 = xpoint[position] - xpoint[position+4];
        x16 = xpoint[position] - xpoint[position+5];


        x23 = xpoint[position+1] - xpoint[position+2];
        x24 = xpoint[position+1] - xpoint[position+3];
        x25 = xpoint[position+1] - xpoint[position+4];
        x26 = xpoint[position+1] - xpoint[position+5];

        x34 = xpoint[position+2] - xpoint[position+3];
        x35 = xpoint[position+2] - xpoint[position+4];
        x36 = xpoint[position+2] - xpoint[position+5];

        x45 = xpoint[position+3] - xpoint[position+4];
        x46 = xpoint[position+3] - xpoint[position+5];

        x56 = xpoint[position+4] - xpoint[position+5];

        double part1, part2, part3, part4, part5, part6;
        part1 = dx2 * dx3 * dx4 * dx5 * dx6 / x12 / x13 / x14 / x15 / x16 * ypoint[position];
        part2 = dx1 * dx3 * dx4 * dx5 * dx6 / (-x12) / x23 / x24 / x25 / x26 * ypoint[position+1];
        part3 = dx1 * dx2 * dx4 * dx5 * dx6 / (-x13) / (-x23) / x34 / x35 / x36 * ypoint[position+2];
        part4 = dx1 * dx2 * dx3 * dx5 * dx6 / (-x14) / (-x24) / (-x34) / x45 / x46 * ypoint[position+3];
        part5 = dx1 * dx2 * dx3 * dx4 * dx6 / (-x15) / (-x25) / (-x35) / (-x45) / x56 * ypoint[position+4];
        part6 = dx1 * dx2 * dx3 * dx4 * dx5 / (-x16) / (-x26) / (-x36) / (-x46) / (-x56) * ypoint[position+5];

        // 	timer::tick("Mathzone","Poly_Inter_xy");
        return part1 + part2 + part3 + part4 + part5 + part6;
    }
    else
    {
        // 	timer::tick("Mathzone","Poly_Inter_xy");
        return ypoint[position];
    }
}

double distance
(
	const Vector3<double> &pos1, 
	const Vector3<double> &pos2, //mohan add 2016-10-21
	const double &a1,
	const double &a2,
	const double &a3
)
{
	const double dx = shortest(pos1.x, pos2.x, a1);
	const double dy = shortest(pos1.y, pos2.y, a2);
	const double dz = shortest(pos1.z, pos2.z, a3);
	return sqrt(dx*dx+dy*dy+dz*dz); 
}

int str2int(const string str)
{
    stringstream ss;
    int y;
    ss<<str;
    ss>>y;
    return y;
}
string int2str(const int y)
{
    stringstream ss;
    string str;
    ss<<y;
    ss>>str;
    return str;
}
string dou2str(const double y)
{
    stringstream ss;
    string str;
    ss<<y;
    ss>>str;
    return str;
}
double str2dou(const string str)
{
    stringstream ss;
    double y;
    ss<<str;
    ss>>y;
    return y;
}
//qianrui 2020-2-18
void searchead(ifstream &ifskwt,string& txt,const string obj,int n)
{
 	int find=0;
	string useless;
	while(ifskwt>>useless)
    {
		getline(ifskwt,txt);
        if(useless==obj)
        {
            find++;
			if(find==n)
				break;
        }
		//cout<<useless<<" "<<txt<<endl;//test
    }
    if(ifskwt.eof())
    {
        cout<<"Can't find  "<<obj<<". File isn't enough!"<<endl;
        exit(0);
    }
	return;
}

double mysecond(long int& time)
{
	return double(time)/CLOCKS_PER_SEC;
}