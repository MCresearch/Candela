#include "HBs.h"
#include "cellFile.h"

#ifndef HB_BREAK_H
#define HB_BREAK_H

class HB_break
{
    public:
    HB_break();
    ~HB_break();

    void Routine();
    void calc(CellFile &cel);
    bool in_list(int &target, int* list, int &list_length);

    double* break_reason_count; // 0 for OO length; 1 for OOH angle; 2 for both;
    double* form_reason_count; // same as above
    int** last_donate_index;
    int** last_donate_H;
    int* last_ndonate;
};

#endif