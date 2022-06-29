#ifndef DATA3D_H
#define DATA3D_H


class Data3D
{
	public: 
	
	Data3D(){};
	~Data3D(){};

	void Routine();
	
	private:

	void Format_PROFESS();
	void Format_clear();
	void Profile(double ***data);

	int nx;
	int ny;
	int nz;

};


#endif
