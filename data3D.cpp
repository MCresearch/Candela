#include "input.h"
#include "data3D.h"

void Data3D::Routine()
{
	if(INPUT.format3D==1)
	{
		Format_clear();
	}
	else if(INPUT.format3D==2)
	{
		Format_PROFESS();
	}
	return;
}

void Data3D::Format_PROFESS()
{
	cout << " Read in the 3D data format : PROFESS " << endl;
	ifstream ifs(INPUT.data_in.c_str());
	if(!ifs)
	{
		cout << " Can't find the file: " << INPUT.data_in << endl;
		exit(0);
	}

	string tmp;
	ifs >> tmp >> this->nx >> tmp >> this->ny >> tmp >> this->nz;
	int nspin;
	ifs >> tmp >> tmp >> tmp >> nspin;

	cout << " dimension of 3D data is " << nx << " " << ny << " " << nz << endl;
	cout << " spin is " << nspin << endl;

	assert(nx>0);
	assert(ny>0);
	assert(nz>0);
	assert(nspin>0);

	// read in the 3D data.
	double*** data;
	data = new double**[nx];
	for(int i=0; i<nx; ++i)
	{
		data[i] = new double*[ny];
		for(int j=0; j<ny; ++j)
		{
			data[i][j] = new double[nz];
		}
	}

	for(int k=0; k<nz; ++k)
	{
		for(int j=0; j<ny; ++j)
		{
	        for(int i=0; i<nx; ++i)
			{
				ifs >> data[i][j][k];
			}
		}
	}

	// output the profile data.
	this->Profile(data);

	// delete the 3D data
	for(int i=0; i<nx; ++i)
	{
		for(int j=0; j<ny; ++j)
		{
			delete[] data[i][j];
		}
		delete[] data[i];
	}
	delete[] data;

	ifs.close();

	return;
}

void Data3D::Format_clear()
{
	// read in the dimension of data.
	ifstream ifs(INPUT.data_in.c_str());
	if(!ifs)
	{
		cout << " Can't find the file:" << INPUT.data_in << endl;
		exit(0);
	}

	string line1, line2, line3;
	getline(ifs, line1);
	getline(ifs, line2);

	string tmp;
	ifs >> tmp >> tmp;

	ifs >> nx >> tmp >> ny >> tmp >> nz;
	getline(ifs, line3);
	
	cout << " dimension is " << nx << " " << ny << " " << nz << endl;

	assert(nx>0);
	assert(ny>0);
	assert(nz>0);


	// read in the 3D data.
	double*** data;
	data = new double**[nx];
	for(int i=0; i<nx; ++i)
	{
		data[i] = new double*[ny];
		for(int j=0; j<ny; ++j)
		{
			data[i][j] = new double[nz];
		}
	}

	for(int k=0; k<nz; ++k)
	{
		for(int j=0; j<ny; ++j)
		{
	        for(int i=0; i<nx; ++i)
			{
				int tmpx, tmpy, tmpz;
				ifs >> tmpx >> tmpy >> tmpz >> data[i][j][k];
			}
		}
	}

	// output the profile data.
	this->Profile(data);


	// delete the 3D data
	for(int i=0; i<nx; ++i)
	{
		for(int j=0; j<ny; ++j)
		{
			delete[] data[i][j];
		}
		delete[] data[i];
	}
	delete[] data;



	return;
}

void Data3D::Profile(double*** data)
{
	double* pdata;
	if(INPUT.direction==1)
	{
		pdata = new double[nx];
	}
	else if(INPUT.direction==2)
	{
		pdata = new double[ny];
	}
	else if(INPUT.direction==3)
	{
		pdata = new double[nz];
	}

	ofstream ofs(INPUT.data_out.c_str());
	if(!ofs)
	{
		cout << " Can't open the file: " << INPUT.data_out << endl;
	}


	if(INPUT.direction==3)
	{
		for(int k=0; k<nz; ++k)
		{
			double v = 0.0;
			for(int i=0; i<nx; ++i)
			{
				for(int j=0; j<ny; ++j)
				{
					v += data[i][j][k];
				}
			}
			ofs << k+1 << " " << v << endl;
		}
	}
	else if(INPUT.direction==2)
	{
		for(int j=0; j<ny; ++j)
		{
			double v = 0.0;
			for(int i=0; i<nx; ++i)
			{
				for(int k=0; k<nz; ++k)
				{
					v += data[i][j][k];
				}
			}
			ofs << j+1 << " " << v << endl;
		}
	}
	else if(INPUT.direction==1)
	{
		for(int i=0; i<nx; ++i)
		{
			double v = 0.0;
			for(int j=0; j<ny; ++j)
			{
				for(int k=0; k<nz; ++k)
				{
					v += data[i][j][k];
				}
			}
			ofs << i+1 << " " << v << endl;
		}
	}
	else
	{
		cout << "Something wrong." << endl;
		exit(0);
	}
	
	ofs.close();

	return;
}
