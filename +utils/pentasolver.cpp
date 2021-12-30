/**
 * @file main.cpp
 * @author theComputeKid
 * @brief Entry point from MATLAB via mex for the penta-diagonal solver.
 * @version 0.1
 * @date 2021-12-28
 *
 * @copyright Copyright (c) 2021 theComputeKid
 */
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include <gpu/mxGPUArray.h>
#include <mex.h>

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cusparse_v2.h>

// #include "mexUtils.hpp"

namespace
{
  struct gpuArray
  {
    void* ptr;
    mxClassID type;
    std::vector<std::size_t> dims;

    gpuArray(mxArray const* const in) : gpuArrayPtr(mxGPUCreateFromMxArray(in))
    {
      type = mxGPUGetClassID(gpuArrayPtr);
      ptr = const_cast<void*>(mxGPUGetDataReadOnly(gpuArrayPtr));
      dims = getDims();
    };

    gpuArray(std::vector<std::size_t> const dims, mxClassID const type)
        : dims(dims), type(type),
          gpuArrayPtr(mxGPUCreateGPUArray(dims.size(), dims.data(), type, mxREAL, MX_GPU_DO_NOT_INITIALIZE))
    {
      ptr = mxGPUGetData(const_cast<mxGPUArray*>(gpuArrayPtr));
    };

    ~gpuArray()
    {
      mxGPUDestroyGPUArray(gpuArrayPtr);
    };

    mxArray* createMxArray()
    {
      return mxGPUCreateMxArrayOnGPU(gpuArrayPtr);
    };

  private:
    mxGPUArray const* const gpuArrayPtr;
    std::vector<size_t> getDims()
    {
      const size_t* dimPtr = mxGPUGetDimensions(gpuArrayPtr);
      const auto nDims = mxGPUGetNumberOfDimensions(gpuArrayPtr);

      if (nDims > 2) mexErrMsgIdAndTxt("pentasolver:invalidInput", "Input array must have a maximum of 2 dimensions.");

      const size_t nRows = *(dimPtr);
      const size_t nCols = *(dimPtr + 1);

      mxFree((void*)dimPtr);

      return {nRows, nCols};
    };
  };

  /**
   * @brief Check if input template is a float.
   *
   * @tparam T Template to check.
   * @return true Is a float.
   * @return false Is not a float.
   */
  template <typename T> bool constexpr iIsFloat()
  {
    return std::is_same_v<T, float>;
  };

  /**
   * @brief Throws an error message and appends the function name to the error ID.
   *
   * @param errID ID is added onto function name with a preceeding semi-colon.
   * @param errMsg Message to display in MATLAB.
   */
  void iThrowErrMsg(std::string const& errID, std::string const& errMsg)
  {
    std::string const funcName = "pentasolver";
    std::string const id = funcName + ":" + errID;
    mexErrMsgIdAndTxt(id.c_str(), errMsg.c_str());
  };

  /**
   * @brief Check the number of inputs are appropriate for the function. Throw error if not.
   *
   * @param actualNumInputs Number of inputs provided.
   * @param expectedNumInputs Number of inputs expected.
   */
  void iCheckNumInputs(std::size_t const actualNumInputs, std::size_t const expectedNumInputs)
  {
    if (actualNumInputs != expectedNumInputs)
    {
      std::string const errorMsg = "Function accepts " + std::to_string(expectedNumInputs) + " inputs but " +
                                   std::to_string(actualNumInputs) + " were given.";
      ::iThrowErrMsg("TooManyInputs", errorMsg.c_str());
    }
  }

  /**
   * @brief Check the number of outputs are not more than appropriate for the function. Throw error if not.
   *
   * @param actualNumOutputs Number of outputs provided.
   * @param expectedNumOutputs Number of outputs requested.
   */
  void iCheckNumOutputs(std::size_t const actualNumOutputs, std::size_t const expectedNumOutputs)
  {
    if (actualNumOutputs > expectedNumOutputs)
    {
      std::string const errorMsg = "Function gives " + std::to_string(expectedNumOutputs) + " outputs but " +
                                   std::to_string(actualNumOutputs) + " were requested.";
      mexErrMsgIdAndTxt("TooManyOutputs", errorMsg.c_str());
    }
  }

  /**
   * @brief Templated helper to get the amount of workspace memory required for the solver.
   *
   * @tparam T
   * @param handle
   * @param dims
   * @param out
   * @param wks
   * @return cusparseStatus_t
   */
  template <typename T>
  void iGetWksEstimate(cusparseHandle_t& handle, std::vector<std::size_t> const& dims, std::vector<T*> const& out,
                       std::size_t& wks)
  {
    std::size_t constexpr algo = 0;
    cusparseStatus_t status;
    if constexpr (::iIsFloat<T>())
    {
      status = cusparseSgpsvInterleavedBatch_bufferSizeExt(handle, algo, dims[0], out[1], out[2], out[3], out[4],
                                                           out[5], out[0], dims[1], &wks);
    }
    else
    {
      status = cusparseDgpsvInterleavedBatch_bufferSizeExt(handle, algo, dims[0], out[1], out[2], out[3], out[4],
                                                           out[5], out[0], dims[1], &wks);
    }
    if (status != CUSPARSE_STATUS_SUCCESS)
      iThrowErrMsg("WksEstimate", "Error in Workspace size estimate: " + std::to_string(status));
  }

  /**
   * @brief Converts a device pointer from aggregate to interleaved format.
   *
   * @tparam T Float or double type.
   * @param in Input device pointer (aggregate form).
   * @param out Output device pointer (interleaved form).
   * @param handle cuBLAS handle.
   * @param N Number of rows of b.
   * @param batchSize Number of columns of b.
   */
  template <typename T>
  void iConvertToInterleaved(std::vector<T*> const& in, std::vector<T*>& out, std::size_t const N,
                             std::size_t const batchSize)
  {
    cublasHandle_t cublasH;
    cublasStatus_t cublasStatus = cublasCreate(&cublasH);

    if (in.size() != out.size())
      iThrowErrMsg("InternalError", "InternalError: number of input and output array pointers do not match.");

    T const one = 1;
    T const zero = 0;

    std::function<cublasStatus_t(cublasHandle_t, cublasOperation_t, cublasOperation_t, int, int, const T*, const T*,
                                 int, const T*, const T*, int, T*, int)>
        xGeamFunc;

    if constexpr (iIsFloat<T>())
    {
      xGeamFunc = cublasSgeam;
    }
    else
    {
      xGeamFunc = cublasDgeam;
    }

    for (std::size_t i = 0; i < in.size(); i++)
    {
      cublasStatus = xGeamFunc(cublasH, CUBLAS_OP_T, /* transa */
                               CUBLAS_OP_T,          /* transb, don't care */
                               batchSize,            /* number of rows of ds */
                               N,                    /* number of columns of ds */
                               &one, in[i],          /* ds0 is n-by-batchSize */
                               N,                    /* leading dimension of ds0 */
                               &zero, NULL, N,       /* don't cae */
                               out[i],               /* ds is batchSize-by-n */
                               batchSize);           /* leading dimension of ds */

      if (cublasStatus != CUBLAS_STATUS_SUCCESS)
        iThrowErrMsg("GEAM", "Error in CUBLAS XGEAM: " + std::to_string(cublasStatus));
    }
    if (cublasH) cublasDestroy(cublasH);
  }

  template <typename T> void iConvertToAggregate(T* const in, T*& out, std::vector<std::size_t> const& dims)
  {
    cublasHandle_t cublasH;
    cublasStatus_t cublasStatus = cublasCreate(&cublasH);

    T const one = 1;
    T const zero = 0;

    std::size_t const N = dims[0];
    std::size_t const batchSize = dims[1];

    std::function<cublasStatus_t(cublasHandle_t, cublasOperation_t, cublasOperation_t, int, int, const T*, const T*,
                                 int, const T*, const T*, int, T*, int)>
        xGeamFunc;

    if constexpr (iIsFloat<T>())
    {
      xGeamFunc = cublasSgeam;
    }
    else
    {
      xGeamFunc = cublasDgeam;
    }

    cublasStatus = xGeamFunc(cublasH, CUBLAS_OP_T, /* transa */
                             CUBLAS_OP_T,          /* transb, don't care */
                             N,                    /* number of columns of ds */
                             batchSize,            /* number of rows of ds */
                             &one, in,             /* ds0 is n-by-batchSize */
                             batchSize,            /* leading dimension of ds0 */
                             &zero, NULL, N,       /* don't cae */
                             out,                  /* ds is batchSize-by-n */
                             N);                   /* leading dimension of ds */

    if (cublasStatus != CUBLAS_STATUS_SUCCESS)
      iThrowErrMsg("GEAM", "Error in CUBLAS XGEAM: " + std::to_string(cublasStatus));
    if (cublasH) cublasDestroy(cublasH);
  }
  /**
   * @brief Templated Solver for GPSV
   *
   * @tparam T
   */
  template <typename T>
  void iExecSolver(cusparseHandle_t& cusparseH, std::vector<std::size_t> const& dims, std::vector<T*> const& out,
                   std::size_t const wks)
  {
    cusparseStatus_t status;
    void* wksPtr;
    cudaError_t const cudaStatus = cudaMalloc((void**)&wksPtr, wks);

    if (cudaStatus != cudaSuccess)
      iThrowErrMsg("WorkspaceAlloc", "Unable to allocate workspace: " + std::to_string(cudaStatus));

    std::size_t constexpr algo = 0;

    if constexpr (iIsFloat<T>())
    {
      status = cusparseSgpsvInterleavedBatch(cusparseH, algo, dims[0], out[1], out[2], out[3], out[4], out[5], out[0],
                                             dims[1], wksPtr);
    }
    else
    {
      status = cusparseDgpsvInterleavedBatch(cusparseH, algo, dims[0], out[1], out[2], out[3], out[4], out[5], out[0],
                                             dims[1], wksPtr);
    }
    cudaFree(wksPtr);
    if (status != CUSPARSE_STATUS_SUCCESS) iThrowErrMsg("GPSV", "Error in GPSV: " + std::to_string(status));
  }

  /**
   * @brief Template specific implementation of the algorithm (Ax=b).
   *
   * @tparam T Template to perform compute for.
   * @param b
   * @param minus2
   * @param minus1
   * @param main
   * @param plus1
   * @param plus2
   * @return ::gpuArray
   */
  template <typename T>
  std::shared_ptr<::gpuArray> iPentasolver(::gpuArray const& b, ::gpuArray const& minus2, ::gpuArray const& minus1,
                                           ::gpuArray const& main, ::gpuArray const& plus1, ::gpuArray const& plus2)
  {

    cusparseHandle_t cusparseH;
    cusparseStatus_t cusparseStatus;
    cusparseStatus = cusparseCreate(&cusparseH);

    // Pointers into which to output interleaved data.
    auto bI = ::gpuArray(b.dims, b.type);
    auto minus2I = ::gpuArray(b.dims, b.type);
    auto minus1I = ::gpuArray(b.dims, b.type);
    auto mainI = ::gpuArray(b.dims, b.type);
    auto plus1I = ::gpuArray(b.dims, b.type);
    auto plus2I = ::gpuArray(b.dims, b.type);

    auto X = std::make_shared<::gpuArray>(b.dims, b.type);

    // Pointers to input diagonal arrays in aggregate format.
    std::vector<T*> in = {(T*)b.ptr, (T*)minus2.ptr, (T*)minus1.ptr, (T*)main.ptr, (T*)plus1.ptr, (T*)plus2.ptr};

    // Pointers to output diagonal arrays in interleaved format.
    std::vector<T*> out = {(T*)bI.ptr, (T*)minus2I.ptr, (T*)minus1I.ptr, (T*)mainI.ptr, (T*)plus1I.ptr, (T*)plus2I.ptr};

    ::iConvertToInterleaved<T>(in, out, b.dims[0], b.dims[1]);

    std::size_t wks = 0;
    ::iGetWksEstimate<T>(cusparseH, b.dims, out, wks);
    ::iExecSolver<T>(cusparseH, b.dims, out, wks);
    T* xPtr = (T*)X->ptr;
    ::iConvertToAggregate<T>((T*)bI.ptr, xPtr, b.dims);
    if (cusparseH) cusparseDestroy(cusparseH);
    return X;
  };
} // namespace

/**
 * @brief MATLAB entry point to penta-diagonal solver based on CUDA.
 */
void mexFunction(int nlhs,        /* number of expected outputs */
                 mxArray* plhs[], /* array of pointers to output arguments */
                 int nrhs,        /* number of inputs */
                 const mxArray* prhs[] /* array of pointers to input arguments */)
{

  ::iCheckNumInputs(nrhs, 6);
  ::iCheckNumOutputs(nlhs, 1);

  auto const b = ::gpuArray(prhs[0]);
  auto const minus2 = ::gpuArray(prhs[1]);
  auto const minus1 = ::gpuArray(prhs[2]);
  auto const main = ::gpuArray(prhs[3]);
  auto const plus1 = ::gpuArray(prhs[4]);
  auto const plus2 = ::gpuArray(prhs[5]);

  if (minus2.type != b.type || minus1.type != b.type || main.type != b.type || plus1.type != b.type ||
      plus2.type != b.type)
    iThrowErrMsg("InconsistentTypes", "Types of all inputs must match.");

  if (b.type != mxSINGLE_CLASS && b.type != mxDOUBLE_CLASS)
    iThrowErrMsg("InvalidDataTypes", "Input type must be single or double.");

  auto x = (b.type == mxSINGLE_CLASS) ? iPentasolver<float>(b, minus2, minus1, main, plus1, plus2)
                                      : iPentasolver<double>(b, minus2, minus1, main, plus1, plus2);

  plhs[0] = x->createMxArray();
}
