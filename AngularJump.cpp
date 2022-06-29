#include "AngularJump.h"
#include "gfun.h"
AngularJump::AngularJump()
{
    this->AJ_ss_index = -1;
    this->AJ_ss_time = -1;
    this->centralO = -1;
    this->Oab = new int[2];
    Oab[0] = -1;
    Oab[1] = -1;
    Oa = -1;
    Ob = -1;
    center_angle = -1;
    H = -1;
}
AngularJump::~AngularJump()
{
    delete[] this->Oab;
}

void AngularJump::setup_AJ(int nAJ, AngularJump* AJ)
{
    ifstream ifs("middle_log_reorg_sort.txt");
    for (int iAJ = 0; iAJ < nAJ; iAJ++)
    {
        ifs >> AJ[iAJ].AJ_ss_index >> AJ[iAJ].AJ_ss_time >> AJ[iAJ].centralO >> AJ[iAJ].Oab[0] >> AJ[iAJ].Oab[1] >> AJ[iAJ].center_angle >> AJ[iAJ].H;
        cout << iAJ << " " << AJ[iAJ].AJ_ss_index << " " << AJ[iAJ].AJ_ss_time << " " << AJ[iAJ].centralO << " " << AJ[iAJ].Oab[0] << " " << AJ[iAJ].Oab[1] << " " << AJ[iAJ].center_angle << " " << AJ[iAJ].H << endl;
    }
    ifs.close();
    return;
}