/*
 * mexTransitions.cc
 *
 *  created: June 2018
 *   author: Kaushik Mallik
 */

#include <iostream>
#include <vector>
/* mex */
#include "mex.h"
#include "ClassHandle.hh"
/* scots */
#include "scots.hh"
// #include "TransitionFunction.hh"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /* get pointer to the command */
  const char* command=mxArrayToString(prhs[0]);
  /* init */
  if(!strcmp(command,"init")) {
    /* get pointer to filename */
    const char* filename=mxArrayToString(prhs[1]);
    /* return a handle to the transition function instance */
    scots::TransitionFunction *tf = new scots::TransitionFunction;
    if(!read_from_file(*tf,filename)) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    plhs[0] = convertPtr2Mat<scots::TransitionFunction>(tf);
    return;
  }

  /* delete */
  if (!strcmp(command,"delete")) {
    /* delete the TransitionFunction */
    destroyObject<scots::TransitionFunction>(prhs[1]);
    return;
  }

  /* get pointer to UniformGrid object */
  scots::TransitionFunction *tf=convertMat2Ptr<scots::TransitionFunction>(prhs[1]);

  /* domain: return states with valid inputs */
  if (!strcmp(command,"domain")) {
    UniformGrid* ss, is;
    // int ab = prhs[1];
    int ab = 0;
    std::string file = "X/X" + std::to_string(ab);
    bool readFlag = read_from_file(*ss, file);
    if (!readFlag)
      error("runtime error: could not read the state space from file");
    else
      int sdim = ss->get_dim();
    bool readFlag = read_from_file(*is, file);
    if (!readFlag)
      error("runtime_error: could not read the input space from file");
    else
      int idim = is->get_dim();

    abs_ptr_type ntr = tf->get_no_transitions();
    double** domain = nullptr;
    domain = new double[ntr][2*sdim+idim];

    tf->get_domain(ss, is, domain);
    /* copy state to std::vector */
    // std::vector<std::vector<double>> domain =
    //   tf->get_domain<std::vector<double>>();
    // /* no. of states */
    // abs_type NS = tf->get_no_states();
    // /* no. of inputs */
    // abs_type NU = tf->get_no_inputs();
    // for (abs_type i=0; i<NS; i++) {
    //
    //   for (abs_type j=0; j<NU; j++) {
    //     std::vector<abs_type> post = tf->get_post();
    //     for (int k=0; k<post.size(); k++)
    //       domain.push_back(post[k]);
    //   }
    // }
    // mwSize N=domain.size();
    if(!ntr) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    // mwSize dim=domain[0].size();
    /* create matrix to store input */
    plhs[0]=mxCreateDoubleMatrix(ntr,(2*sdim+idim),mxREAL);
    double **dom=mxGetPr(plhs[0]);
    for(abs_ptr_type i=0; i<ntr; i++) {
      for(int j=0; j<2*sdim+idim; j++) {
        dom[i][j]=domain[i][j];
      }
    }
    return;
  }
  return;
}
