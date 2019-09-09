#ifndef BLACKBOXREACHWRAPPER_HH_
#define BLACKBOXREACHWRAPPER_HH_
/*
 * BlackBoxReachWrapper.hh
 *
 *  created on: 06.09.2019
 *      author: kaushik
 */

/*
 * A wrapper function to set up the task of black box abstr refinement
 *
 */

#include <array>
#include <iostream>
#include <cmath>
#include <time.h> /* time used to seed the random number */
#include <float.h> /* for the smallest positive number */
#define _USE_MATH_DEFINES

#include "BlackBoxReach.hh"

using namespace std;
using namespace scots;
using namespace helper;

/*! Convert a vector of arrays to a vector of vectors of the same size */
template<std::size_t SIZE, class T>
void vecArr2vecVec(const std::vector<std::array<T,SIZE>> va, std::vector<std::vector<T>>& vv) {
    vv.clear();
    for (int l=0; l<va.size(); l++) {
        std::vector<T> d(std::begin(va[l]), std::end(va[l]));
        vv.push_back(d);
    }
}

/*! Simulate the system using the controller already available in the BlackBoxReach object */
template<std::size_t SIZE_o, std::size_t SIZE_g, std::size_t SIZE_i, class sys_type, class X_type, class U_type>
void simulateSystem(BlackBoxReach* abs,
                    const std::vector<std::array<double,SIZE_o>> ho,
                    const std::vector<std::array<double,SIZE_g>> hg,
                    const std::vector<std::array<double,SIZE_i>> hi,
                    sys_type sys_post,
                    X_type x,
                    U_type u,
                    std::vector<double>& unsafeAt,
                    std::vector<std::vector<double>>& sys_traj) {
//    double toss1, toss2;
    std::vector<std::vector<double>> hovec, hgvec;
//    /* pick a random initial point within the provided initial set */
//    toss1 = 0.01*(rand() % (int)(100*(hi[0][1]+hi[0][0])))-hi[0][0];
//    toss2 = 0.01*(rand() % (int)(100*(hi[0][3]+hi[0][2])))-hi[0][2];
//    std::vector<double> init = {toss1, toss2};
//    sys_traj.push_back(init);
    vecArr2vecVec(ho,hovec);
    vecArr2vecVec(hg,hgvec);
    abs->simulateSys(sys_post, x, u, sys_traj,hovec,hgvec,unsafeAt);
}

/*! Compute the size of default array */
template<class T>
inline int size_of_array(T* arr) {
    return (sizeof(arr)/sizeof(arr[0]));
}

/*! Compute the multi-layered abstraction of the systemm. */
template<class X_type, class U_type, class sys_type, class rad_type, class O_type, class G_type, class I_type, class HO_type, class HG_type, class HI_type, class ho_type, class hg_type, class hi_type>
double find_abst(X_type x, U_type u,
              sys_type sys_post, rad_type radius_post,
              double* lbX, double* ubX, double* lbU, double* ubU,
              double* etaX, double* etaU, double tau, double systemTau,
              int numAbs, double* etaRatio, double tauRatio,
              O_type spawnO, G_type spawnG, I_type spawnI,
              HO_type HO, HG_type HG, HI_type HI,
              ho_type ho, hg_type hg, hi_type hi,
              int nSubInt, int systemNSubInt, int p,
              int NN, X_type explRadius, double reqd_success_rate,
              bool readTsFromFile, bool useColors, const char* logfile, int verbose=0) {
    
    int dimX = x.size();
    int dimU = u.size();
    
    int act_success_count = 0;
    bool reqd_success_rate_reached = false;
    bool spec_with_no_refinement;
    
    System sys(dimX, lbX, ubX, etaX, tau,
                   dimU, lbU, ubU, etaU,
                   etaRatio, tauRatio, nSubInt, numAbs);
    
    BlackBoxReach* abs= new BlackBoxReach(logfile,verbose);
    abs->initialize(&sys,systemTau,systemNSubInt);
    if (readTsFromFile) {
        abs->loadTs();
    } else {
        /* start with only coarse transitions */
        abs->initializeAbstraction(sys_post,radius_post,x,u);
        checkMakeDir("T");
        saveVec(abs->Ts_, "T/T");
        clog << "Wrote Ts_ to file.\n";
    }
    
    /* start with a fresh copy of the abstraction */
//    BlackBoxReach* abs = new BlackBoxReach(*abs_ref);
    /* we keep a moving maximum */
    double spec = 0;
    /* spec_old is the older value computed: at first it is set to -1. */
    double spec_old = -1;
    /* initialize variables */
    double toss1, toss2, toss3;
    std::vector<double> unsafeAt;
    std::vector<std::vector<double>> sys_traj, abs_traj;
    double distance;
    int iter=1;
    while (1) {
        cout << "\033[1;4;34m\n\nIteration = "<< iter <<" \n\n\033[0m";
        /*************************************************************/
        /* **** Computation of SPEC with abstraction refinement **** */
        /*************************************************************/
        cout << "\033[1;4;34mStarting computation of SPEC.\n\n\033[0m";
        spec_with_no_refinement=true;
        /* iterate over all the environments */
        for (int e=0; e<NN; e++) {
            /* spawn environment */
            if (!useColors)
                cout << "Environment #" << e << "\n";
            spawnO(HO,ho,0.0,verbose);
            spawnG(HG,hg,0.0,verbose);
            spawnI(HI,hi,verbose);
            
            /* *** synthesize controller on the available abstraction only *** */
            /* initialize the environment in the abstraction */
            bool flag = abs->initializeSpec(HO,ho,HG,hg,HI,hi);
            if (!flag) { /* ignore this environment */
                abs->clear();
                continue;
            }
            if (verbose>0)
                cout << "Abstraction initialized with the specification.\n";
            /* do a simple multi-layered reachability (without further refinement) */
            abs->plainReach(p);
            if (verbose>0)
                cout << "Controller synthesis done.\n";
            /* if the abstract game is not winning, the distance remains 0 */
            distance = 0;
            /* check if there is a controller for the abstraction */
            if (abs->isInitWinning()) {
                if (verbose>0)
                    cout << "There is a controller.\n";
                /* simulate the concrete system using the synthesized controller */
                // debug
//                abs->writeVecToFile(sys_traj,"Figures/traj.txt");
                /* pick a random initial point within the provided initial set */
                double toss1, toss2;
                toss1 = 0.01*(rand() % (int)(100*(hi[0][1]+hi[0][0])))-hi[0][0];
                toss2 = 0.01*(rand() % (int)(100*(hi[0][3]+hi[0][2])))-hi[0][2];
                std::vector<double> init = {toss1, toss2};
                sys_traj.push_back(init);
                simulateSystem(abs,ho,hg,hi,sys_post,x,u,unsafeAt,sys_traj);
                //debug
//                abs->writeVecToFile(sys_traj,"Figures/sys_traj.txt","clean");
//                abs->saveFinalResult();
                //debug end
//                printArray(sys_traj, )
                /* check if the system trajectory was unsafe */
                if (unsafeAt.size()!=0) {
                    if (useColors) {
                        /* print in red (the code 31) */
                        cout << "\033[31mEnvironment #" << e << "\033[0m\n";
                    }
                    if (verbose>0) {
                        cout << "System trajectory went to unsafe part at (";
                        for (int i=0; i<dimX; i++) {
                            cout << unsafeAt[i] << ",";
                        }
                        cout << "\b).\n";
                    }
                    /* simulate the abstract trajectory and measure the minimum distance from the obstacles */
                    abs_traj.push_back(sys_traj[0]);
                    std::vector<std::vector<double>> hovec;
                    vecArr2vecVec(ho,hovec);
                    abs->simulateAbs(abs_traj,hovec,distance);
                    cout << "Distance of abstract trajectory from the unsafe states = " << distance << ".\n";
                    //debug
//                    abs->writeVecToFile(abs_traj,"Figures/abs_traj.txt","clean");
                    //debug end
                    if (distance>spec) {
                        /* Refine abstraction if unsafe until the specified accuracy is reached */
                        //debug
//                        checkMakeDir("T");
//                        saveVec(abs_ref->Ts_, "T/T");
                        //debug end
                        while (1) {
                            cout << "\tStarting a refinement loop to minimize the distance.\n";
                            abs->clear();
                            bool flag1 = abs->exploreAroundPoint(unsafeAt, explRadius, u, sys_post, radius_post);
                            if (!flag1) { /* already explored upto the finest layer */
                                cout << "\tNo more refinement possible. Final distance = " << distance << "\n";
                                break;
                            }
                            // debug
//                            checkMakeDir("T");
//                            saveVec(abs->Ts_, "T/T");
                            //end
                            /* even if the abstract game solving was good enough with the previous spec value, still then the abstract games need to be solved again, as the abstraction went finer during the spec computation */
                            spec_with_no_refinement = false;
                            
                            /* initialize distance; distance remains 0 if the abstract game with the refined abstraction fails */
                            distance = 0;
                            /* do a simple multi-layered reachability (without further refinement) */
                            abs->initializeSpec(HO,ho,HG,hg,HI,hi);
                            abs->plainReach(p);
                            //debug
//                            abs->saveFinalResult();
                            //debug end
                            if (!abs->isInitWinning()) {
                                cout << "\tUpdated distance = 0\n";
                                break;
                            }
                            unsafeAt.clear();
                            /* simulate the system using the synthesized controller */
                            /* use the old (problematic) initial state */
                            std::vector<double> old_init = sys_traj[0];
                            sys_traj.clear();
                            sys_traj.push_back(old_init);
                            simulateSystem(abs,ho,hg,hi,sys_post,x,u,unsafeAt,sys_traj);
                            //debug
//                            abs->writeVecToFile(sys_traj,"Figures/sys_traj.txt","clean");
                            //debug end
                            if (unsafeAt.size()==0) {
                                if (verbose>0)
                                    cout << "\tThe controller worked for the system as well.\n";
                                cout << "\tUpdated distance = 0\n";
                                break;
                            }
                            /* if the controller failed on the system, recompute distance and continue refining the abstraction around the point of failure, which is stored in the variable unsafeAt */
                            abs_traj.clear();
                            abs_traj.push_back(sys_traj[0]);
                            abs->simulateAbs(abs_traj,hovec,distance);
                            cout << "\tUpdated distance = " << distance << "\n";
                            /* if the distance after refinement became smaller than the current spec value, no need to refine further */
                            if (distance<=spec) {
                                break;
                            }
                            //debug
//                            abs->writeVecToFile(abs_traj,"Figures/abs_traj.txt","clean");
                            //debug end
                            //debug
//                            checkMakeDir("T");
//                            saveVec(abs_ref->Ts_, "T/T");
//                            abs->saveFinalResult();
                            //debug end
                            
                        } /* End of while loop */
                    } /* End of if (distance > spec) */
                    
                    /* update SPEC */
                    if (distance>spec) {
                        spec = distance;
                    }
                } else { /* if (unsafeAt.size()==0) */
                    if (useColors) {
                        /* print in green (the code 32) */
                        cout << "\033[32mEnvironment #" << e << "\033[0m\n";
                    }
                    if (verbose>0)
                        cout << "The controller worked for the system as well.\n";
                } /* End of if (unsafeAt.size()==0) */
            } else { /* if the abstract initial states are not winning */
                if (useColors)
                    cout << "Environment #" << e << "\n";
                if (verbose>0)
                    cout << "There is no controller for this environment.\n";
            } /* End of if(abs->isInitWinning()) */
            
//            if (distance>0)
//                act_success_count[0]++;
            
            /* clear the last environment */
            ho.clear();
            hg.clear();
            hi.clear();
            HO.clear();
            HG.clear();
            HI.clear();
            /* reset the controller synthesis related members of the abstraction object */
            abs->clear();
            /* clear the sets used for trajectory simulation */
            unsafeAt.clear();
            sys_traj.clear();
            abs_traj.clear();
        } /* End of for loop over the set of environments */
        cout << "\n\nSPEC = " << spec << "\n\n";
        
        /* if the last computed spec value worked well enough in terms of abstract game solving, and the spec value didn't increase from the last computed spec value, and finally if the abstraction was not refined further during the computation of spec, then we are done */
        if ((reqd_success_rate_reached) &&
            (spec <= spec_old) &&
            (spec_with_no_refinement)) {
            cout << "Abstract computation finished.\nTermination condition used: reqd. abstract synthesis success rate reached + the SPEC value didn't increase afterwards + no refinement was performed during the last SPEC loop.\n";
            /* write outputs to files */
            checkMakeDir("T");
            saveVec(abs->Ts_, "T/T");
            cout << "\nWrote transitions to files. Exiting process.\n\n";
            break;
        }
        
        /************************************************************************/
        /* Abstract game solving with refinement: Refinement of abstraction to maximize number of winning environments */
        /************************************************************************/
        cout << "\033[1;4;34mStarting abstraction refinement for environment satisfaction.\n\n\033[0m";
        /* keep a history of the 'seen' environments */
        std::vector<std::array<SymbolicSet*,3>> ENV_HIST;
        /* initial abstraction */
//        BlackBoxReach* abs = new BlackBoxReach(*abs_ref);
        int unique_env_count = 0; /* number of environments explored excluding duplicate cases */
        act_success_count = 0; /* reset the actual success count for solving the abstract game */
        /* how far the environments are to be made conservative (no role for computing SPEC) */
        double eps = DBL_MIN; /* very small number added to make sure that boundary cases are pessimistically resolved */
        double spec2 = spec + eps;
        /* iterate over the environments */
        for (int e=0; e<NN; e++) {
            abs->clear();
            /* spawn environment */
            if (!useColors)
                cout << "Environment #" << e << "\n";
            spawnO(HO,ho,spec2,verbose);
            spawnG(HG,hg,spec2,verbose);
            spawnI(HI,hi,verbose);
            
            /* initialize the environment in the abstraction */
            bool flag = abs->initializeSpec(HO,ho,HG,hg,HI,hi);
            if (!flag) /*ignore this specificaiton */
                continue;
            if (verbose>0)
                cout << "Abstraction initialized with the specification.\n";
            /* if this is not the fist environment, then assume that the environment is seen before */
            bool newenv;
            if (e==0) {
                newenv = true;
            } else {
                newenv = false;
            }
            for (int i=0; i<ENV_HIST.size(); i++) {
                /* check similarity in the obstacles in the finest layer */
                if (abs->Os_[numAbs-1]->symbolicSet_!=ENV_HIST[i][0]->symbolicSet_) {
                    newenv = true;
                }
                /* check similarity in the goals in the finest layer */
                if (abs->Gs_[numAbs-1]->symbolicSet_!=ENV_HIST[i][1]->symbolicSet_) {
                    newenv = true;
                }
                /* check similarity in the initial states in the finest layer */
                if (abs->X0s_[numAbs-1]->symbolicSet_!=ENV_HIST[i][2]->symbolicSet_) {
                    newenv = true;
                }
                
                if (!newenv) /* similarity found */
                    break;
                else if (i<ENV_HIST.size()-1) /* if this is not the last iteration of this for loop then reset newenv flag and continue with the next one */
                    newenv = false;
            }
            /* if the envrionment is seen before, continue with the next one */
            if (!newenv) {
                if (verbose>0)
                    cout << "Environment seen before. Coninuing with the next one.\n";
                continue;
            }
            unique_env_count++;
            
            /* store the environment info in the list */
            SymbolicSet* O = new SymbolicSet(*abs->Os_[numAbs-1]);
            SymbolicSet* G = new SymbolicSet(*abs->Gs_[numAbs-1]);
            SymbolicSet* X0 = new SymbolicSet(*abs->X0s_[numAbs-1]);
            std::array<SymbolicSet*,3> arr = { O, G, X0};
            ENV_HIST.push_back(arr);
            
            /* perform a lazy reach-avoid abstraction refinement step */
            abs->onTheFlyReach(p, sys_post, radius_post, x, u);
            /* print result of the abstraction refinement */
            if (abs->isInitWinning()) {
                /* increment the success rate counter */
                act_success_count++;
                if (useColors)
                /* print in green (the code 32) */
                    cout << "\033[32mEnvironment #" << e << "\033[0m\n";
                if (verbose>0)
                    cout << "The environment is now winning.\n";
            } else {
                if (useColors)
                /* print in red (the code 31) */
                    cout << "\033[31mEnvironment #" << e << "\033[0m\n";
                if (verbose>0)
                    cout << "The environment was not winnable.\n";
            }
            
            /* clear the specification sets */
            ho.clear();
            hg.clear();
            hi.clear();
            HO.clear();
            HG.clear();
            HI.clear();
        } /* End of for loop over the set of environments */
        if (unique_env_count==0) {
            cout << "\nSPEC value "<< spec <<" is too high. The multi-layered abstraction is too coarse. Consider recomputation with 1 more layer at the bottom.";
            return (-1);
        }
        if ((double)act_success_count/unique_env_count > reqd_success_rate) {
            cout << "\nTarget precision reached.";
//            break;
            /* the computation is finished, provided the SPEC value doesn't increase */
            reqd_success_rate_reached = true;
        }
        if (abs->Ts_[numAbs-1]->symbolicSet_==abs->ddmgr_->bddOne()) {
            cout << "\nNo more refinement possible.";
            /* write outputs to files */
            checkMakeDir("T");
            saveVec(abs->Ts_, "T/T");
            cout << "\nWrote transitions to files. Exiting process.\n\n";
            break;
        }
        /* clear the controller related vaiables */
        abs->clear();
//        BlackBoxReach* abs_ref = new BlackBoxReach(*abs);
        /* store the current spec value and reset the old spec variable */
        spec_old = spec;
        spec = 0;
        
        iter++;
//        for (int l=0; l<dimX; l++)
//        act_success_count = 0;
    } /* End of the big while loop over spec computation and game solving */
    return spec;
}

/*! Test the multi-layered abstraction with a number of randomly generated test environments from a given set of environments */
template<class X_type, class U_type, class sys_type, class rad_type, class O_type, class G_type, class I_type, class HO_type, class HG_type, class HI_type, class ho_type, class hg_type, class hi_type>
void test_abstraction(BlackBoxReach* abs, double spec_final,
                      X_type x, U_type u,
                      sys_type sys_post, rad_type radius_post,
                      O_type spawnO, G_type spawnG, I_type spawnI,
                      HO_type HO, HG_type HG, HI_type HI,
                      ho_type ho, hg_type hg, hi_type hi,
                      int p,
                      int NN, int num_tests,
                      bool useColors, const char* logfile, int verbose=0) {
    cout << "\033[1;4;34mStarting to test the computed abstraction for controller synthesis against test environments.\n\n\033[0m";
    int dimX = x.size();
    /* perform a series of tests on the computed abstract transition system */
    int success_count = 0;
    int num_cont_found = 0;
    for (int e=0; e<num_tests; e++) {
        bool validEnv = false;
        while (!validEnv) {
            spawnO(HO,ho,spec_final,verbose);
            spawnG(HG,hg,spec_final,verbose);
            spawnI(HI,hi,verbose);
            validEnv = abs->initializeSpec(HO,ho,HG,hg,HI,hi);
        }
        if (verbose>0)
            cout << "Abstraction initialized with the specification.\n";
        /* do a simple multi-layered reachability (without further refinement) */
        abs->plainReach(p);
        // debug
//        abs->onTheFlyReach(p,sys_post,radius_post,x,u);
//        abs->saveFinalResult();
        // debug end
        if (verbose>0)
            cout << "Controller synthesis done.\n";
        if (abs->isInitWinning()) {
            if (verbose>0)
                cout << "There is a winning controller.\n";
            clog << "\nFinal number of controllers: " << abs->finalCs_.size() << '\n';
            /* now try the controller on the system itself */
            std::vector<std::vector<double>> sys_traj;
            /* pick a random initial point within the provided initial set */
            double toss1, toss2;
            toss1 = 0.01*(rand() % (int)(100*(hi[0][1]+hi[0][0])))-hi[0][0];
            toss2 = 0.01*(rand() % (int)(100*(hi[0][3]+hi[0][2])))-hi[0][2];
            std::vector<double> init = {toss1, toss2};
            sys_traj.push_back(init);
            std::vector<double> unsafeAt;
            simulateSystem(abs,ho,hg,hi,sys_post,x,u,unsafeAt,sys_traj);
            // debug
//            abs->writeVecToFile(sys_traj,"Figures/sys_traj.txt","clean");
//            abs->saveFinalResult();
            // debug end
            if (unsafeAt.size()!=0) {
                if (useColors) {
                    /* print in red (the code 31) */
                    cout << "\033[31mEnvironment #" << e << "\033[0m\n";
                }
                if (verbose>0) {
                    cout << "System trajectory went to unsafe part at (";
                    for (int i=0; i<dimX; i++) {
                        cout << unsafeAt[i] << ",";
                    }
                    cout << "\b).\n";
                }
            } else {
                if (useColors) {
                    /* print in green (the code 32) */
                    cout << "\033[32mEnvironment #" << e << "\033[0m\n";
                }
                if (verbose>0)
                    cout << "The controller worked for the system as well.\n";
                success_count++;
            }
        } else {
            if (verbose>0)
                cout << "No controller found for this environment.\n";
            cout << "Environment #" << e << "\n";
            num_cont_found++;
        }
        /* clear the specification sets */
        ho.clear();
        hg.clear();
        hi.clear();
        HO.clear();
        HG.clear();
        HI.clear();
        /* clear the controller info */
        abs->clear();
        
        //        /* save synthesis results */
        //        abs->saveFinalResult();
    }
    cout << "\nTest result:\n";
    if (num_tests!=num_cont_found) {
        cout << "Success rate = " << success_count*100/(num_tests-num_cont_found) << "%.\n";
        if (num_cont_found!=0) {
            cout << "This excludes the " << num_cont_found*100/num_tests << "% cases where no controller could be synthesized.";
        }
    }
    else
        cout << "All tests failed.\n";
}

#endif
