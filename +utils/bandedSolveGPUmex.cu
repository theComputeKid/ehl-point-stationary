#include <array>
#include <cstddef>
#include <string>

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse.h>

#include "mex.h"

#include "gpu/mxGPUArray.h"

namespace
{
  void iMakeInterleaved(void const* const inPtr, void*& outPtr, std::size_t const batchSize,
                        std::size_t const numBatches, mxClassID const type)
  {
    std::size_t const elemBytes = (type == mxSINGLE_CLASS) ? sizeof(float) : sizeof(double);
    auto const cudaStatus = cudaMalloc((void**)&outPtr, elemBytes * batchSize * numBatches);

    if (cudaStatus != cudaSuccess)
    {
      char constexpr errId[] = "parallel:gpu:bandedSolveGPUmex:MallocError";
      char constexpr errMsg[] = "Unable to allocate GPU memory";
      mexErrMsgIdAndTxt(errId, errMsg);
    }

    float constexpr one = 1;
    float constexpr zero = 0;
    cublasHandle_t cublasH = NULL;

    auto cublasStatus = cublasCreate(&cublasH);
    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
    {
      char constexpr errId[] = "parallel:gpu:bandedSolveGPUmex:InitError";
      std::string const errMsg = "Unable to init cublas: Error ";
      mexErrMsgIdAndTxt(errId, (errMsg + std::to_string(cublasStatus)).c_str());
    }

    cublasStatus = cublasSgeam(cublasH, CUBLAS_OP_T,    /* transa */
                               CUBLAS_OP_T,             /* transb, don't care */
                               batchSize,               /* number of rows of ds */
                               numBatches,              /* number of columns of ds */
                               &one, (float*)inPtr,     /* ds0 is n-by-batchSize */
                               numBatches,              /* leading dimension of ds0 */
                               &zero, NULL, numBatches, /* don't cae */
                               (float*)outPtr,          /* ds is batchSize-by-n */
                               batchSize);              /* leading dimension of ds */

    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
    {
      char constexpr errId[] = "parallel:gpu:bandedSolveGPUmex:InterleaveError";
      std::string const errMsg = "Unable to interleave (cublasXgeam): Error ";
      mexErrMsgIdAndTxt(errId, (errMsg + std::to_string(cublasStatus)).c_str());
    }
  };

} // namespace

// Penta-Diagonal Banded Solver
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, mxArray const* prhs[])
{
  std::size_t constexpr N = 5;          // Number of diagonals
  std::size_t constexpr inArgs = N + 1; // Number of inputs

  // Check number of inputs
  if (nrhs != inArgs)
  {
    char constexpr errId[] = "parallel:gpu:bandedSolveGPUmex:InvalidInputs";
    char constexpr errMsg[] = "Invalid number of inputs.";
    mexErrMsgIdAndTxt(errId, errMsg);
  }

  // Ensure everything is a GPUarray
  for (auto i = 0; i < nrhs; i++)
    if (!mxIsGPUArray(prhs[i]))
    {
      char constexpr errId[] = "parallel:gpu:bandedSolveGPUmex:NotOnGPU";
      std::string const errMsg = "Input not a GPU array: ";
      mexErrMsgIdAndTxt(errId, (errMsg + std::to_string(i)).c_str());
    }

  /* Initialize the MathWorks GPU API. */
  mxInitGPU();

  // Extract GPUArrays from all inputs
  std::array<mxGPUArray const* const, 5> const bands = {mxGPUCreateFromMxArray(prhs[0]),  // minus 2
                                                        mxGPUCreateFromMxArray(prhs[1]),  // minus 1
                                                        mxGPUCreateFromMxArray(prhs[2]),  // main
                                                        mxGPUCreateFromMxArray(prhs[3]),  // plus 1
                                                        mxGPUCreateFromMxArray(prhs[4])}; // plus 2
  auto const Y = mxGPUCreateFromMxArray(prhs[5]);                                         // RHS

  // Ensure everything is the same data type
  auto const dataType = mxGPUGetClassID(Y);
  for (auto i = 0; i < N; i++)
  {
    if (mxGPUGetClassID(bands[i]) != dataType)
    {
      char constexpr errId[] = "parallel:gpu:bandedSolveGPUmex:InconsistentType";
      char constexpr errMsg[] = "Inputs are of different types";
      mexErrMsgIdAndTxt(errId, errMsg);
    }
  }

  // Get raw pointers
  std::array<void const*, 5> const bandsPtr = {mxGPUGetDataReadOnly(bands[0]),  // minus 2
                                               mxGPUGetDataReadOnly(bands[1]),  // minus 1
                                               mxGPUGetDataReadOnly(bands[2]),  // main
                                               mxGPUGetDataReadOnly(bands[3]),  // plus 1
                                               mxGPUGetDataReadOnly(bands[4])}; // plus 2
  void const* const YPtr = mxGPUGetDataReadOnly(Y);                             // RHS

  // Step 1: Convert to Interleave
  std::array<void*, 5> bandsInterleavedPtr;
  void* YInterleavedPtr;

  // Total elements
  auto const dims = mxGPUGetDimensions(Y);
  auto const batchesSize = dims[0];
  auto const numBatches = dims[1];
  mxFree((void*)dims);

  for (auto i = 0; i < N; i++)
  {
    iMakeInterleaved(bandsPtr[i], bandsInterleavedPtr[i], batchesSize, numBatches, dataType);
    mxGPUDestroyGPUArray(bands[i]);
  }

  iMakeInterleaved(YPtr, YInterleavedPtr, batchesSize, numBatches, dataType);
  mxGPUDestroyGPUArray(Y);

  for (auto i = 0; i < N; i++)
  {
    cudaFree(bandsInterleavedPtr[i]);
  }
  cudaFree(YInterleavedPtr);

  //   for (auto i = 0; i < nrhs)
  //     auto const ds = ;
  //   auto const dl = mxGPUCreateFromMxArray(prhs[0]);
  //   auto const ds = mxGPUCreateFromMxArray(prhs[0]);
  //   auto const ds = mxGPUCreateFromMxArray(prhs[0]);
  //   auto const ds = mxGPUCreateFromMxArray(prhs[0]);
  //   auto const ds = mxGPUCreateFromMxArray(prhs[0]);

  //   /* Declare all variables.*/
  //   mxGPUArray const* A;
  //   mxGPUArray* B;
  //   double const* d_A;
  //   double* d_B;
  //   int N;

  //   /* Choose a reasonably sized number of threads for the block. */
  //   int const threadsPerBlock = 256;
  //   int blocksPerGrid;
}
