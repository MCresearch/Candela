#include <complex>
using namespace std;

void multipAHB(int M, int N, int K, 
    complex<double> *A, int LDA,
    complex<double> *B, int LDB,
    complex<double> *C, int LDC);

void dtrimultipAHB(int M, int N, int K, 
    complex<double> *A, int LDA,
    complex<double> *B, int LDB,
    complex<double> *C, int BIAS);