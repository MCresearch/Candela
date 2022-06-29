#include "input.h"
#include "pseudo.h"

void Pseudo::Routine()
{
	TITLE("Pseudo","Routine");

	cal();

	return;
}


void Pseudo::cal()
{
	TITLE("Pseudo","cal");
	
	const double pi = 3.1415926535897;

	ifstream ifs(INPUT.pseudo_in.c_str());
	if(!ifs)
	{
		cout << " Can not find the file: " << INPUT.pseudo_in << endl;
		exit(0);
	}
	else
	{
		cout << " Find the file: " << INPUT.pseudo_in << endl;
	}

	// get first line;
	string line;
	std::getline(ifs,line);
	std::getline(ifs,line);
	std::getline(ifs,line);
	std::getline(ifs,line);
	cout << " Ignore the first 4 lines" << endl;

	// get second line;
	double v1,v2;
	ifs >> v1 >> v2; 
	std::getline(ifs,line);

	// get third line;
	double Gmax;
	ifs >> Gmax;
	cout << " Gmax = " << Gmax << endl;
	assert(Gmax>0.0);

	// get the local pseudopotentials unitl the mark '1000'
	int nline=0;
	double a,b,c;
	while(!ifs.eof() )
	{
		ifs >> a; 
		if( abs(a-1000.0) < 1.0e-5 ) break;
		ifs >> b;
		READ_VALUE(ifs, c);
		++nline;
	}
	cout << " Nubmer of lines : " << nline << endl;




	//---------------------------------------
	// read in the pseudopotentials. 
	//---------------------------------------
	if(!ifs.good())
	{
		cout << " Something is wrong during reading" << endl;
		cout << " Check the last number of pseudopotential" << endl;
		exit(0);
	}
	cout << " Rewind the readin pointer." << endl;
	string line2;
	assert(nline>0);
	// to the top again.
	ifs.seekg(0);
	std::getline(ifs,line2);
	std::getline(ifs,line2);
	std::getline(ifs,line2);
	std::getline(ifs,line2);
	double t1,t2,t3;
	ifs >> t1 >> t2 >> t3;
	cout << " trash=" << t1 << " " << t2 << " " << t3 << endl;


	double* pp = new double [3*nline];
	for(int i=0; i<3*nline; ++i)
	{
		ifs >> pp[i];
	}
	ifs.close();
	
	// output the pseudopotential file.
	ofstream ofs(INPUT.pseudo_out.c_str());
	if(!ofs)
	{
		cout << " Can not open file : " << INPUT.pseudo_out << endl;
		exit(0);
	}
	else
	{
		assert(Gmax > 0.0);
		double dg = Gmax / (nline*3-1) * 0.529177; // mohan fix bug 2013-03-30
		double mev_ang3 = 13.6058*2*0.529177*0.529177*0.529177;
		int z = INPUT.pseudo_z;
		cout << " z = " << z << endl;
		if(INPUT.pseudo_z<=0)
		{
			cout << " check the input value of pseudo_z" << endl;
			exit(0);
		}
		cout << " the output data format : G, vloc(G) [contains -4*PI*Z/G^2], vloc(G)+4*PI*Z/G^2 [pure local LPS]" << endl;
		for(int i=1; i<3*nline; ++i)
		{
			double g_norm=dg*(double)i;
			ofs << g_norm << " " << pp[i]/mev_ang3 << " " << pp[i]/mev_ang3 + 4.0*pi*(double)z/g_norm/g_norm << endl;
		}
	}
	ofs.close();

	delete[] pp;
}
