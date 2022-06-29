#include "DoubleDonor.h"
#include "gfun.h"
DoubleDonor::DoubleDonor()
{
    this->DD_start_time = -1;
    this->DD_end_time = -1;
    this->centralO = -1;
    this->Oab = new int[2];
    Oab[0] = -1;
    Oab[1] = -1;
    Obefore = -1;
    Oafter = -1;
    H = -1;
    jump = "False";
}
DoubleDonor::~DoubleDonor()
{
    delete[] this->Oab;
}

void DoubleDonor::setup_DD(int &nDD, DoubleDonor* DD)
{
    ifstream ifs("double_donor.txt");
    string useless;
    for (int i=0; i<10; i++)
    {
        ifs >> useless;
    }
    int iDD2=0;
    double start_time;
    double end_time;
    int centralO;
    int H;
    int* Oab = new int[2];
    int Obefore;
    int Oafter;
    string jump;
    for (int iDD1 = 0; iDD1 < nDD; iDD1++)
    {
        ifs >> start_time >> end_time >> centralO >> H >> Oab[0] >> Oab[1] >> Obefore >> Oafter >> jump >> useless;
        if (iDD1%INPUT.nbin == 0 and Obefore>=0 and Oafter >=0)
        {
            DD[iDD2].DD_start_time = start_time;
            DD[iDD2].DD_end_time = end_time;
            DD[iDD2].centralO = centralO;
            DD[iDD2].H = H;
            DD[iDD2].Oab[0] = Oab[0];
            DD[iDD2].Oab[1] = Oab[1];
            DD[iDD2].Obefore = Obefore;
            DD[iDD2].Oafter = Oafter;
            DD[iDD2].jump = jump;
            cout << iDD2 << " " << DD[iDD2].DD_start_time << " " << DD[iDD2].DD_end_time << " " << DD[iDD2].centralO << " " << DD[iDD2].H << " " << DD[iDD2].Oab[0] << " " << DD[iDD2].Oab[1] << " " << DD[iDD2].Obefore << " " << DD[iDD2].Oafter << " " << DD[iDD2].jump << endl;
            iDD2++;
        }
    }
    nDD = iDD2;
    ifs.close();
    //cout << "setup done." << endl;
    return;
}