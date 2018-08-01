/*
 * mexSymbolicSet.cc
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


void check_arguments(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  check_arguments(nlhs, plhs, nrhs, prhs);
  /* get pointer to the command */
  const char* command=mxArrayToString(prhs[0]);
  /* init */
  if(!strcmp(command,"init")) {
    /* get pointer to filename */
    const char* filename=mxArrayToString(prhs[1]);
    /* initialize cudd manager and load BDD and SymbolicSet from file */ 
    Cudd* mgr = new Cudd;
    BDD* bdd = new BDD;
    scots::SymbolicSet *set = new scots::SymbolicSet;
    if(!scots::read_from_file(*mgr,*set,*bdd,filename)) {
      mexErrMsgTxt("mexSymbolicSet: could not read SymbolicSet from filename.");
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      plhs[1]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    /* store ptr into return value */
    plhs[0]=mxCreateNumericMatrix(3,1, mxUINT64_CLASS, mxREAL);
    uint64_t *p=(uint64_t*)mxGetData(plhs[0]);
    p[0]=convertPtr2MatUINT64<scots::SymbolicSet>(set);
    p[1]=convertPtr2MatUINT64<BDD>(bdd);
    p[2]=convertPtr2MatUINT64<Cudd>(mgr);
    /* dimension */
    plhs[1]=mxCreateDoubleMatrix(1,1, mxREAL);
    double *d=mxGetPr(plhs[1]);
    d[0]=set->get_dim();
    return;
  }

  /* check pointer to object handels */
  uint64_t* p=(uint64_t*) mxGetData(prhs[1]);

  /* delete */
  if (!strcmp(command,"delete")) {
    /* delete the objects */
    destroyObject<scots::SymbolicSet>(p[0]);
    destroyObject<BDD>(p[1]);
    destroyObject<Cudd>(p[2]);
    return;
  }

  /* extract pointer information */
  scots::SymbolicSet *set = convertMat2Ptr<scots::SymbolicSet>(p[0]);
  BDD* bdd = convertMat2Ptr<BDD>(p[1]);
  Cudd* mgr = convertMat2Ptr<Cudd>(p[2]);

  /* return grid points */ 
  if(!strcmp(command,"gridpoints")) {
    /* get projection indices */
    double* indices=mxGetPr(prhs[2]);
    const mwSize *dims=mxGetDimensions(prhs[2]);
    std::vector<int> project_dim {}; 
    for(int i=0; i<dims[0]; i++) 
      project_dim.push_back(indices[i]);

    auto grid_points = set->projection(*mgr,*bdd,project_dim);
    size_t no_grid_points = grid_points.size()/dims[0];
    /* convert abstract cell IDs to the center of cells */
    plhs[0]=mxCreateDoubleMatrix(no_grid_points,dims[0],mxREAL);
    double *ptr=mxGetPr(plhs[0]);
    for(size_t i=0; i<no_grid_points; i++) {
      for(int j=0; j<dims[0]; j++) {
        ptr[j*no_grid_points+i]=grid_points[i*dims[0]+j];
      }
    }
    return;
  } 

  /* return grid points restricted to x in rhs[1] */ 
  if(!strcmp(command,"restriction")) {
    /* get restriction state */
    double* xp=mxGetPr(prhs[2]);
    const mwSize *dims=mxGetDimensions(prhs[2]);
    double* ip=mxGetPr(prhs[3]);
    std::vector<double> xv {}; 
    std::vector<int> iv {}; 
    for(int i=0; i<dims[0]; i++) {
      xv.push_back(xp[i]);
      iv.push_back(ip[i]);
    }
    size_t codim = set->get_dim()-dims[0];
    auto grid_points = set->restriction(*mgr,*bdd,xv,iv);
    plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
    if(!grid_points.size() || !codim) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    size_t no_grid_points = grid_points.size()/codim;
    /* convert abstract cell IDs to the center of cells */
    plhs[0]=mxCreateDoubleMatrix(no_grid_points,codim,mxREAL);
    double *ptr=mxGetPr(plhs[0]);
    for(size_t i=0; i<no_grid_points; i++) {
      for(int j=0; j<codim; j++) {
        ptr[j*no_grid_points+i]=grid_points[i*codim+j];
      }
    }
    return;
  } 
  return;
}


void check_arguments(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs < 1) 
    mexErrMsgTxt("mexSymbolicSet: arguments missing.");
  /* input must be a string */
  if (mxIsChar(prhs[0]) != 1)
    mexErrMsgTxt("mexSymbolicSet: first parameter must be a string out of {init,delete,gridpoints,restriction}.");
  const char* command=mxArrayToString(prhs[0]);
  if(!strcmp(command,"init")) {
    /* second argument must also be a character */
    if (nrhs < 2) 
      mexErrMsgTxt("mexSymbolicSet('init',filename): filename missing");
    if (mxIsChar(prhs[1]) != 1)
      mexErrMsgTxt("mexSymbolicSet('init',filename): filename must be a string containing the filename");
    return;
  }
  if(!strcmp(command,"delete")) {
    if (nrhs < 2) 
      mexErrMsgTxt("mexSymbolicSet('delete',obj.h): obj_h missing");
    if (mxGetClassID(prhs[1])!=mxUINT64_CLASS)
      mexErrMsgTxt("mexSymbolicSet('delete',obj.h): second input must the obj.h (a uint64_t array).");
    const mwSize *dims=mxGetDimensions(prhs[1]);
    if (dims[0]!=3 || dims[1]!=1)
      mexErrMsgTxt("mexSymbolicSet('delete',obj.h): size(obj.h) not equal to [3 1].");
      
  }
  if(!strcmp(command,"gridpoints")) {
    if (nrhs < 3) 
      mexErrMsgTxt("mexSymbolicSet('gridpoints',obj.h,p_dim): one of the arguments missing");
    if (mxGetClassID(prhs[1])!=mxUINT64_CLASS)
      mexErrMsgTxt("mexSymbolicSet('gridpoints',obj.h,p_dim): second input must the obj.h (a uint64_t array).");
    const mwSize *dims=mxGetDimensions(prhs[1]);
    if (dims[0]!=3 || dims[1]!=1)
      mexErrMsgTxt("mexSymbolicSet('gridpoints',obj.h,p_dim): size(obj.h) not equal to [3 1].");
  }
  if(!strcmp(command,"restriction")) {
    if (nrhs < 4) 
      mexErrMsgTxt("mexSymbolicSet('restriction',obj.h,x,r_dim): one of the arguments missing");
    if (mxGetClassID(prhs[1])!=mxUINT64_CLASS)
      mexErrMsgTxt("mexSymbolicSet('restriction',obj.h,x,r_dim): second input must the obj.h (a uint64_t array).");
    const mwSize *dims=mxGetDimensions(prhs[1]);
    if (dims[0]!=3 || dims[1]!=1)
      mexErrMsgTxt("mexSymbolicSet('restriction',obj.h,x,r_dim): size(obj.h) not equal to [3 1].");

    const mwSize *dx=mxGetDimensions(prhs[2]);
    const mwSize *di=mxGetDimensions(prhs[3]);
    if(dx[0]!=di[0] || dx[1]!=di[1])
      mexErrMsgTxt("mexSymbolicSet('restriction',obj.h,x,r_dim): size(x)!=size(r_dim).");
  }
}
