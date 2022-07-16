#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <cassert>
#include <ctime>

using namespace std;

int main()
{
	//------------------------------------------------------
	//
	//  PLEASE SETUP BY YOURSELF FOR DIFFERENT SITUATION
	//
	// prepare data:
	// how many processos are used (split_to_nfile)
	int split_to_nfile=10;
	// how many configurations (max)
	int max=50000;
	// how manly lines (line)
	int line=115;
	// which directory we want to begin with
	int begin_from_dir=2;

	//------------------------------------------------------

	// static structure fafctor (ssf)
	// reciprocal space vector (g)
	double *ssf = new double[line];
	double *g = new double[line];
	for(int i=0; i<line; ++i) ssf[i] = 0.0;
	for(int i=0; i<line; ++i) g[i] = 0.0;
	double tmp_d, tmp_g;

	// reading geometry
	int count=0;
	for(int dir=begin_from_dir; dir<split_to_nfile; ++dir)
	{
		stringstream ss;
		ss << "cal_ssf_" << dir << "/Li_ssf.txt";
		ifstream ifs(ss.str().c_str());

		if(ifs)
		{
			cout << " Read data from " << ss.str() << endl;
			for(int i=0; i<line; ++i)
			{
				ifs >> tmp_g >> tmp_d;
				g[i] = tmp_g;
				ssf[i] += tmp_d;
			}

			ifs.close();
			++count;
		}
		else
		{
			cout << " can't find the file." << endl;
			exit(0);
		}
	}

	// output data
	//
	// get averaged static structure factor
	cout << " count files found =" << count << endl;
	assert(count>0);
	for(int i=0; i<line; ++i)
	{
		ssf[i] /= count;
	}

	string output_file = "Final_ssf.txt";
	ofstream ofs(output_file.c_str());
	for(int i=0; i<line; ++i)
	{
		ofs << g[i] << " " << ssf[i] << endl;
	}
	ofs.close();
	cout << " Final data is in " << output_file << endl;

	return 0;
}
