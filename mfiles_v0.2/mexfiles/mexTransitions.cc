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
    // std::cout << "debug 1" << '\n';
    scots::UniformGrid ss, is;
    // UniformGrid* is;
    int ab = (int)*mxGetPr(prhs[2]);
    // int ab = 0;
    int sdim, idim;
    std::string file = "X/X" + std::to_string(ab-1);
    // const char* filename=mxArrayToString(prhs[2]);
    // std::cout << "debug 2" << '\n';
    bool readFlag = read_from_file(ss, file.c_str());
    // std::cout << "debug 2" << '\n';
    if (!readFlag)
      error("runtime error: could not read the state space from file");
    // // else
    sdim = ss.get_dim();
    readFlag = read_from_file(is, "U");
    // std::cout << "debug 2" << '\n';
    if (!readFlag)
      error("runtime_error: could not read the input space from file");
    // // else
    idim = is.get_dim();

    abs_ptr_type ntr = tf->get_no_transitions();
    // std::cout << "debug 2" << '\n';
    if(!ntr) {
      plhs[0]=mxCreateDoubleMatrix(0,0,mxREAL);
      return;
    }
    double** domain = new double*[ntr];//[2*(mwSize)sdim+(mwSize)idim];
    for (size_t i = 0; i < ntr; i++) {
      domain[i] = new double[2*sdim+idim];
    }
    // std::cout << "debug 3" << '\n';
    tf->get_domain(&ss, &is, domain);
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

    // std::cout << "debug 4" << '\n';
    // mwSize dim=domain[0].size();
    /* create matrix to store input */
    plhs[0]=mxCreateDoubleMatrix(ntr,(2*sdim+idim),mxREAL);
    // std::cout << "debug 5" << '\n';
    double *dom=mxGetPr(plhs[0]);
    // std::cout << "debug 6" << '\n';
    for(abs_ptr_type i=0; i<ntr; i++) {
      // std::cout << "debug 7 \t " << i << '\n';
      for(int j=0; j<2*sdim+idim; j++) {
        // std::cout << "debug 8 \t " << j  << '\n';
        dom[j*ntr+i]=domain[i][j];
      }
    }
    return;
  }
  return;
}
