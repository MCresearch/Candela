#include "matrixmultip.h"
#include <iostream>
#include <chrono>

static double multi_cost = 0;
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
#ifdef USE_CUDA
void gpu_dtrimultipahb(int m, int n, int k,
                       complex<double> *a, int lda,
                       complex<double> *b, int ldb,
                       complex<double> *c, int bias);
#endif

void multipAHB(int M, int N, int K,
               complex<double> *A, int LDA,
               complex<double> *B, int LDB,
               complex<double> *C, int LDC)
{
    auto start_time = std::chrono::high_resolution_clock::now();
    multipahb_(M,N,K,A,LDA,B,LDB,C,LDC);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    multi_cost += (double)(elapsed.count()) / 1000.0;
    return;
}

void dtrimultipAHB(int M, int N, int K,
                   complex<double> *A, int LDA,
                   complex<double> *B, int LDB,
                   complex<double> *C, int BIAS,
                   const string &device)
{
    auto start_time = std::chrono::high_resolution_clock::now();
#ifdef USE_CUDA
    if (device == "auto" || device == "gpu")
    {
        gpu_dtrimultipahb(M, N, K, A, LDA, B, LDB, C, BIAS);
    }
    else
#endif
        dtrimultipahb_(M, N, K, A, LDA, B, LDB, C, BIAS);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    multi_cost += (double)(elapsed.count()) / 1000.0;

    return;
}

void multi_time()
{
    std::cout << "Matrix multiplication time elapsed: " << multi_cost << "s" << std::endl;
}