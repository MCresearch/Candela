#ifndef STRESS_H
#define STRESS_H



class stress_average
{
public:
	stress_average();
	~stress_average();

	void Routine();

	double* stress_list;
};

#endif