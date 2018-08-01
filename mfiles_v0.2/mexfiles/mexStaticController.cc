/*
 * mexStaticController.cc
 *
 *  created: Jan 2017
 *   author: Matthias Rungger
 */

#include <iostream>
#include <vector>
/* mex */
#include "mex.h"
#include "ClassHandle.hh"
/* scots */
#include "scots.hh"



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* get pointer to the command */
  const char* command=mxArrayToString(prhs[0]);
  /* init */
  if(!strcmp(command,"init")) {
    /* get pointer to filename */
    const char* filename=mxArrayToString(prhs[1]);
    /* return a handle to the uniform grid instance */
    scots::StaticController *controller = new scots::StaticController;
    if(!scots::read_from_file(*controller,filename)) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    plhs[0] = convertPtr2Mat<scots::StaticController>(controller);
    return;
  }
  /* delete */
  if (!strcmp(command,"delete")) {
    /* delete the StaticController */
    destroyObject<scots::StaticController>(prhs[1]);
    return;
  }

  /* get pointer to UniformGrid object */ 
  scots::StaticController *controller=convertMat2Ptr<scots::StaticController>(prhs[1]);

  /* control: return control input associated with x */
  if (!strcmp(command,"control")) {
    /* copy state to std::vector */
    double *px=mxGetPr(prhs[2]);
    const mwSize *dims=mxGetDimensions(prhs[2]);
    std::vector<double> x(dims[0]);
    for(size_t i=0; i<dims[0]; i++) {
      x[i]=px[i];
    }
    /* get input associated with x */
    try {
      std::vector<std::vector<double>> 
      u=controller->get_control<std::vector<double>,std::vector<double>>(x);
      /* create matrix to store input */
      plhs[0]=mxCreateDoubleMatrix((mwSize)u.size(),(mwSize)u[0].size(),mxREAL);
      double *inputs=mxGetPr(plhs[0]);
      for(size_t i=0; i<u.size(); i++) {
        for(size_t j=0; j<u[i].size(); j++) {
          inputs[j*u.size()+i]=u[i][j];
        }
      }
      return;
    } catch(...) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
  }

  /* domain: return states with valid inputs */
  if (!strcmp(command,"domain")) {
    /* copy state to std::vector */
    std::vector<std::vector<double>> domain =
      controller->get_domain<std::vector<double>>();
    mwSize N=domain.size();
    if(!N) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    mwSize dim=domain[0].size();
    /* create matrix to store input */
    plhs[0]=mxCreateDoubleMatrix(N,dim,mxREAL);
    double *dom=mxGetPr(plhs[0]);
    for(size_t i=0; i<N; i++) {
      for(size_t j=0; j<dim; j++) {
        dom[j*N+i]=domain[i][j];
      }
    }
    return;
  }
  return;
}

