#include "matrixmultip.h"

extern "C"{
    void multipahb_(int& M, int& N, int& K, 
    complex<double> *A, int& LDA,
    complex<double> *B, int& LDB,
    complex<double> *C, int& LDC);

    void dtrimultipahb_(int& M, int& N, int& K, 
    complex<double> *A, int& LDA,
    complex<double> *B, int& LDB,
    complex<double> *C, int& BIAS);
}

void multipAHB(int M, int N, int K, 
    complex<double> *A, int LDA,
    complex<double> *B, int LDB,
    complex<double> *C, int LDC)
{
    multipahb_(M,N,K,A,LDA,B,LDB,C,LDC);
    return;
}

void dtrimultipAHB(int M, int N, int K, 
    complex<double> *A, int LDA,
    complex<double> *B, int LDB,
    complex<double> *C, int BIAS)
{
    dtrimultipahb_(M,N,K,A,LDA,B,LDB,C,BIAS);
    return;
}
