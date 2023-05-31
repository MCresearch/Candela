#ifndef __PWTEST
#define __PWTEST
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include <iostream>
using namespace std;
namespace UT
{
    extern int NPROC, RANK;
}

class TestEnv : public testing::Environment
{
public:
    virtual void SetUp()
    {
        if (UT::RANK == 0)
        {
            cout << "\033[32m"
                 << "[ SET UP UNIT TESTS OF CANDELA]"
                 << "\033[0m" << endl;
            cout << "\033[32m[ " << UT::NPROC << " processors are used."
                 << " ]\033[0m" << endl;
        }
    }
    virtual void TearDown()
    {
        if (UT::RANK == 0)
        {
            cout << "\033[32m"
                 << "[ TEAR DOWN TESTS ]"
                 << "\033[0m" << endl;
        }
    }
};

#endif