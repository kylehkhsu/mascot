/*
 * mexGoal.cc
 *
 *  created: July 2018
 *   author: Kaushik Mallik
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
  // int h = mexPrintf("Entered the function \n");
  // h = mexPrintf("Everything good.\n");
  const char* command=mxArrayToString(prhs[0]);
  /* init */
  if(!strcmp(command,"init")) {/* get pointer to filename */
    const char* filename=mxArrayToString(prhs[1]);
    /* return a handle to the uniform grid instance */
    scots::Goal *target = new scots::Goal;
    if(!scots::read_from_file(*target,filename)) {
      // h = mexPrintf("could not read from the file");
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      // plhs[0]=mxCreateDoubleScalar(10); // debug purpose
      return;
      // h = mexPrintf("not returned");
    }
    plhs[0] = convertPtr2Mat<scots::Goal>(target);
    return;
  }

  /* delete */
  if (!strcmp(command,"delete")) {
    /* delete the StaticController */
    destroyObject<scots::Goal>(prhs[1]);
    return;
  }

  /* get pointer to Goal object */
  scots::Goal *target=convertMat2Ptr<scots::Goal>(prhs[1]);

  /* parameters */
  if (!std::strcmp(command,"parameters")) {
    /* dim */
    size_t dim = target->state_grid.get_dim();
    if(!dim) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    plhs[0]=mxCreateDoubleScalar((double)dim);
    /* eta */
    if(nlhs>=1) {
      std::vector<double> v = target->state_grid.get_eta();
      plhs[1]=mxCreateDoubleMatrix(1,(mwSize)dim,mxREAL);
      double *eta=mxGetPr(plhs[1]);
      for (size_t i = 0; i < v.size(); i++) {
        eta[i] = v[i];
      }
    }
  }

  /* domain: return states which are in the goal */
  if (!strcmp(command,"points")) {
    // h=  mexPrintf("Entered points.\n");
    /* copy state to std::vector */
    std::vector<std::vector<double>> domain =
      target->get_domain<std::vector<double>>();
    mwSize N=domain.size();
    if(!N) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    mwSize dim=domain[0].size();
    /* create matrix to store the states */
    plhs[0]=mxCreateDoubleMatrix(N,dim,mxREAL);
    double *dom=mxGetPr(plhs[0]);
    for(size_t i=0; i<N; i++) {
      for(size_t j=0; j<dim; j++) {
        dom[j*N+i]=domain[i][j];
      }
    }
    // h = mexPrintf("Everything good.\n");
    // nlhs=1;
    return;
  }
  /* checkstate: checks whether a given state is in the goal or not */
  if (!strcmp(command,"checkstate")) {
    mwSize nx    = mxGetM(prhs[2]);
    double* ptr_x = mxGetPr(prhs[2]);
    plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
    double *inside = mxGetPr(plhs[0]);

    std::vector<double> x;
    // std::vector<size_t> idx;
    for(size_t i=0; i<nx; i++) {
      x.push_back(ptr_x[i]);
      // idx.push_back(i);
    }
    if (target->check_state(x)) {
      int h = mexPrintf("is inside\n");
      inside[0] = 1;
    }
    else {
      int h = mexPrintf("is not inside\n");
      inside[0] = 0;
    }
    return;
  }
  return;
}
