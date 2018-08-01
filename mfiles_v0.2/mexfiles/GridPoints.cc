/*
 * GridPoints.cc
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */
#include <iostream>
#include <vector>
/* mex */
#include "mex.h"
/* scots */
#include "scots.hh"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* there must be one argument and it must be a string */
  if (nrhs != 1 || mxIsChar(prhs[0]) != 1) {
    mexErrMsgIdAndTxt("MATLAB:GridPoints",
                      "requires at least one input argument and it must be string.");
  }
  /* get pointer to the filename */
  const char* filename=mxArrayToString(prhs[0]);

  /* read grid info from file */
  scots::UniformGrid grid;
  if(!scots::read_from_file(grid,filename)) {
    mexErrMsgIdAndTxt("MATLAB:GridPoints",
                      "couldn't not read UniformGrid from file: %s ",filename);
  }

  /* read grid point indices */
	bool success=true;
  scots::FileReader reader(filename);
  if(!reader.open()) {
    success=false;
  }
  std::vector<scots::abs_type> abs_gp {};
  if(!reader.get_VECTOR(SCOTS_GP_DATA,abs_gp)) {
    success=false;
  }
	reader.close();
  if(!success) {
    mexErrMsgIdAndTxt("MATLAB:GridPoints",
                      "couldn't not read grid points from file: %s",filename);
	}
  mwSize dim = grid.get_dim();
  mwSize no_grid_points = abs_gp.size();

  /* convert abstract cell IDs to the center of cells */
  plhs[0]=mxCreateDoubleMatrix(no_grid_points,dim,mxREAL);
  double *ptr=mxGetPr(plhs[0]);
  std::vector<double> x(grid.get_dim());
  for(size_t i=0; i<no_grid_points; i++) {
    grid.itox(abs_gp[i],x);
    for(size_t j=0; j<dim ; j++) {
      ptr[j*no_grid_points+i]=x[j];
    }
  }
  return;
}

