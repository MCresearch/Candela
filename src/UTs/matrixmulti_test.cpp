#include "unittest.h"
#include "../matrixmultip.h"
#include <complex>
#include <cstdlib>

TEST(matrixmultip, dtrimultipAHB)
{
    constexpr int M = 40;
    constexpr int N = 40;
    constexpr int K = 450;
    constexpr int LDA = 500;
    constexpr int LDB = 500;
    constexpr int BIAS = 4;

    std::complex<double> a[LDA * M];
    std::complex<double> b[LDB * N];
    for (int i = 0; i < LDA * M; ++i)
    {
        double r = std::rand() / double(RAND_MAX);
        double arg = std::rand() / double(RAND_MAX) * 2 * M_PI;
        a[i] = r * std::complex<double>{cos(arg), sin(arg)};
    }
    for (int i = 0; i < LDB * N; ++i)
    {
        double r = std::rand() / double(RAND_MAX);
        double arg = std::rand() / double(RAND_MAX) * 2 * M_PI;
        b[i] = r * std::complex<double>{cos(arg), sin(arg)};
    }
    constexpr int DimC = M - BIAS - N + 1 > 0 ? (M - BIAS) * N - N * (N - 1) / 2 : (M - BIAS + 1) * (M - BIAS) / 2;
    std::complex<double> c[DimC];
    std::complex<double> expected_c[DimC];

    for (int j = 0, q = 0; j < N; ++j)
    {
        for (int i = j + BIAS; i < M; ++i, ++q)
        {
            std::complex<double> temp = {0.0, 0.0};
            for (int l = 0; l < K; ++l)
            {
                temp += conj(a[i * LDA + l]) * b[j * LDB + l];
            }
            expected_c[q] = temp;
        }
    }

    // Call the cpu function
    dtrimultipAHB(M, N, K, a, LDA, b, LDB, c, BIAS, "cpu");

    // Check the result
    for (int i = 0; i < DimC; i++)
    {
        EXPECT_NEAR(c[i].real(), expected_c[i].real(), 1e-8);
        EXPECT_NEAR(c[i].imag(), expected_c[i].imag(), 1e-8);
    }
#ifdef USE_CUDA
    std::complex<double> c_gpu[DimC];
    // Call the gpu function
    dtrimultipAHB(M, N, K, a, LDA, b, LDB, c_gpu, BIAS, "gpu");
    // Check the result
    for (int i = 0; i < DimC; i++)
    {
        EXPECT_NEAR(c_gpu[i].real(), expected_c[i].real(), 1e-8);
        EXPECT_NEAR(c_gpu[i].imag(), expected_c[i].imag(), 1e-8);
    }
#endif

    multi_time();
}