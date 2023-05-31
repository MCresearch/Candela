#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <complex>
#include <iostream>

__global__ void dtrimultipahb(int m, int n, int k, cuDoubleComplex *a, int lda, cuDoubleComplex *b, int ldb, cuDoubleComplex *c, int bias)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= m || j >= n || i < j + bias)
    {
        return;
    }

    cuDoubleComplex temp = make_cuDoubleComplex(0.0, 0.0);

    for (int l = 0; l < k; l++)
    {
        temp = cuCadd(temp, cuCmul(cuConj(a[i * lda + l]), b[j * ldb + l]));
    }

    int ij = j * (m - bias - 1) - (j - 1) * j / 2 + i - bias;
    c[ij] = temp;
}

void gpu_dtrimultipahb(int m, int n, int k, std::complex<double> *a, int lda, std::complex<double> *b, int ldb, std::complex<double> *c, int bias)
{

    cuDoubleComplex *d_a, *d_b, *d_c;
    cudaMalloc(&d_a, sizeof(cuDoubleComplex) * lda * m);
    cudaMalloc(&d_b, sizeof(cuDoubleComplex) * ldb * n);
    int dimc = m - bias - n + 1 > 0 ? (m - bias) * n - n * (n - 1) / 2 : (m - bias + 1) * (m - bias) / 2;
    cudaMalloc(&d_c, sizeof(cuDoubleComplex) * dimc);

    cudaMemcpy(d_a, a, sizeof(cuDoubleComplex) * lda * m, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, sizeof(cuDoubleComplex) * ldb * n, cudaMemcpyHostToDevice);

    dim3 block_size(16, 16);
    dim3 num_blocks((m + block_size.x - 1) / block_size.x, (n + block_size.y - 1) / block_size.y);

    dtrimultipahb<<<num_blocks, block_size>>>(m, n, k, d_a, lda, d_b, ldb, d_c, bias);

    cudaMemcpy(c, d_c, sizeof(cuDoubleComplex) * dimc, cudaMemcpyDeviceToHost);

    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
}
