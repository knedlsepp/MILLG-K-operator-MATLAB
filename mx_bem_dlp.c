/**
 * @brief Contains MATLAB bindings to the bem_matlab module to compute
 *        the double-layer potential for Laplace in three dimensions.
 */

#include "mex.h"
#include "bem_matlab.h"

#define ERROR_BUFFER_LENGTH 1024

void mexFunction(int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs)
{
  const char* functionName = mexFunctionName();
  char errorMessage[ERROR_BUFFER_LENGTH];

  int nC = 0, nE = 0;
  double* dlp = NULL;
  const double* coordinates = NULL;
  const double* elements    = NULL;
  const double* areas       = NULL;
  double minus_id_factor    = 0.;

  if (nlhs > 1 || (nrhs != 3 && nrhs != 4))
  {
    sprintf( errorMessage,
              "Invalid number of input/output arguments in %s.\n"
              "Please call this function using:\n\n"
              "  K = %s(coordinates,elements,areas[,minus_id_factor=0])\n\n"
              "where coordinates and elements describe the mesh, areas "
              "contains the\n"
              "areas of the respective elements and minus_id_factor is a "
              "factor.\n\n"
              "Then this function compute (K - minus_id_factor), where K\n"
              "denotes the double layer potential for Laplace in three "
              "dimensions.\n\n\n"
			  "Seppcomment: \n"
			  "All elements have to be oriented OUTWARDS and then we use: \n"
			  "bdOp = bsxfun(@rdivide,-buildK(bdcoordinates,boundary,areas,0.5),areas)\n"
			  "for the strayfield.\n"
			  "This may be due to buildK implementing ( (K-1/2), \\chi_j  )_{L^2}\n"
			  "or something similar.",
              functionName, functionName );
    mexErrMsgTxt( errorMessage );
  }

  coordinates = (const double*) mxGetPr(prhs[0]);
  elements    = (const double*) mxGetPr(prhs[1]);
  areas       = (const double*) mxGetPr(prhs[2]);

  if (nrhs == 4)
  {
    minus_id_factor = *((const double*) mxGetPr(prhs[3]));
  }

  nC = mxGetM(prhs[0]);
  nE = mxGetM(prhs[1]);

  plhs[0] = mxCreateDoubleMatrix(nE, nC, mxREAL);
  dlp     = mxGetPr(plhs[0]);

  build_dlp_for_C(
    nC,
    coordinates,
    nE,
    elements,
    areas,
    minus_id_factor,
    dlp);
}

