#ifndef BLACKBOXREACHWRAPPER_HH_
#define BLACKBOXREACHWRAPPER_HH_
/*
 * BlackBoxReachWrapper.hh
 *
 *  created on: 06.09.2019
 *      author: kaushik
 */

/*
 * A set of wrapper functions to set up the task of black box abstr refinement
 *
 */

#include <array>
#include <iostream>
#include <cmath>
#include <time.h> /* time used to seed the random number */
#include <float.h> /* for the smallest positive number */
#include <omp.h> /* for paralle implementation */
#define _USE_MATH_DEFINES

#define CHUNK_SIZE 4

#include "BlackBoxReach.hh"

//using namespace std;
using namespace scots;
using namespace helper;

struct closed_loop_log {
    std::vector<int> abstraction_used;
    std::vector<std::vector<double>> trajectory;
    std::vector<std::vector<double>> strategy;
};

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
template<std::size_t SIZE_o, std::size_t SIZE_g, std::size_t SIZE_i, class sys_type, class X_type, class U_type, class L>
void simulateSystem(BlackBoxReach* abs,
                    const std::vector<std::array<double,SIZE_o>> ho,
                    const std::vector<std::array<double,SIZE_g>> hg,
                    const std::vector<std::array<double,SIZE_i>> hi,
                    sys_type sys_post,
                    X_type x,
                    U_type u,
                    std::vector<double>& unsafeAt,
                    L& sys_log) {
//    double toss1, toss2;
    std::vector<std::vector<double>> hovec, hgvec;
//    /* pick a random initial point within the provided initial set */
//    toss1 = 0.01*(rand() % (int)(100*(hi[0][1]+hi[0][0])))-hi[0][0];
//    toss2 = 0.01*(rand() % (int)(100*(hi[0][3]+hi[0][2])))-hi[0][2];
//    std::vector<double> init = {toss1, toss2};
//    sys_traj.push_back(init);
    vecArr2vecVec(ho,hovec);
    vecArr2vecVec(hg,hgvec);
    abs->simulateSys(sys_post, x, u, sys_log,hovec,hgvec,unsafeAt);
}

/*! Compute the size of default array */
template<class T>
inline int size_of_array(T* arr) {
    return (sizeof(arr)/sizeof(arr[0]));
}

/*! Compute the multi-layered abstraction of the systemm. */
template<class X_type, class U_type, class sys_type, class rad_type , class O_type, class G_type, class I_type, class HO_type, class HG_type, class HI_type, class ho_type, class hg_type, class hi_type, class gen_init_type>
double find_abst(X_type x, U_type u,
              sys_type sys_post, rad_type radius_post,
              double* lbX, double* ubX, double* lbU, double* ubU,
              double* etaX, double* etaU, double tau, double systemTau,
              int numAbs, double* etaRatio, double tauRatio,
              O_type spawnO, G_type spawnG, I_type spawnI,
              HO_type HO, HG_type HG, HI_type HI,
              ho_type ho, hg_type hg, hi_type hi, gen_init_type generateInitial,
              const int nSubInt, const int systemNSubInt, const int p,
              const int NN, const X_type explRadius, const double explHorizon, const double* reqd_success_rate, const double spec_max,
              bool readTsFromFile, bool useColors, const char* logfile, const int verbose=0) {
    
    int dimX = x.size();
    int dimU = u.size();
    
    int act_success_count[2] = {0, 0};
    bool reqd_success_rate_reached = false;
//    bool spec_with_no_refinement;
    
    System sys(dimX, lbX, ubX, etaX, tau,
                   dimU, lbU, ubU, etaU,
                   etaRatio, tauRatio, nSubInt, numAbs);
    
    BlackBoxReach* abs= new BlackBoxReach(logfile,verbose);
    abs->initialize(&sys,systemTau,systemNSubInt);
    // debug
//    abs->Xs_[1]->printInfo(1);
    // debug end
    if (readTsFromFile) {
        abs->loadInitTs();
    } else {
        /* start with only coarse transitions along with the auxiliary exploration abstractions */
        abs->initializeAbstractionWithExplore(sys_post,radius_post,x,u);
        checkMakeDir("T_init");
        saveVec(abs->Ts_, "T_init/T");
        clog << "Wrote Ts_ to file.\n";
        checkMakeDir("uT_init");
        saveVec(abs->uTs_, "uT_init/uT");
        clog << "Wrote uTs_ to file.\n";
    }
    
    /* start with a fresh copy of the abstraction */
//    BlackBoxReach* abs = new BlackBoxReach(*abs_ref);
    /* we keep a moving maximum */
    double spec = 0;
    /* spec_old is the older value computed: at first it is set to same as spec. */
    double spec_old = 0;
    /* initialize variables */
    double toss1, toss2, toss3;
    std::vector<double> unsafeAt;
    
    closed_loop_log sys_log, abs_log;
    double distance;
    int iter=1;
    while (1) {
        cout << "\033[1;4;34m\n\nIteration = "<< iter <<" \n\n\033[0m";
        /************************************************************************/
        /* Abstract game solving with refinement: Refinement of abstraction to maximize number of winning environments */
        /************************************************************************/
        cout << "\033[1;4;34mStarting abstraction refinement for environment satisfaction.\n\n\033[0m";
        /* keep a history of the 'seen' environments */
        std::vector<std::array<SymbolicSet*,3>> ENV_HIST;
        /* initial abstraction */
//        BlackBoxReach* abs = new BlackBoxReach(*abs_ref);
        int unique_env_count = 0; /* number of environments explored excluding duplicate cases */
        act_success_count[1] = 0; /* reset the actual success count for solving the abstract game */
        /* how far the environments are to be made conservative (no role for computing SPEC) */
//        double eps = 1e-13; /* very small number added to make sure that boundary cases are pessimistically resolved */
        double spec2 = spec + eps;

        /* create all NN number of separate BlackBoxReach objects with the same variables as abs */
//        BlackBoxReach** abs_copies= new BlackBoxReach*[NN];
        std::vector<BlackBoxReach*> abs_copies;
        for (int e=0; e<NN; e++) {
            BlackBoxReach* abs_copy = new BlackBoxReach(*abs);
//            abs_copies[e] = abs_copy;
            abs_copies.push_back(abs_copy);
        }
#pragma omp parallel
        {
#pragma omp parallel for \
schedule(dynamic, CHUNK_SIZE) \
shared(abs,ENV_HIST,spec2,unique_env_count) private(HO,HG,HI,ho,hg,hi)
        /* iterate over the environments */
        for (int e=0; e<NN; e++) {
            BlackBoxReach* abs_curr;
#pragma omp critical (get_abs_copy)
            {
                abs_curr = abs_copies.back();
                abs_copies.pop_back();
            }
            if (verbose>=0) {
                std::cout << "Starting computation for environment: " << e << ".\n";
            }
            /* spawn environment */
            if (!useColors)
                cout << "Environment #" << e << "\n";
            spawnO(HO,ho,verbose);
            spawnG(HG,hg,verbose);
            spawnI(HI,hi,verbose);
            
            /* initialize the environment in the abstraction */
            bool flag;
//#pragma omp critical (initialization)
            flag = abs_curr->initializeSpec(HO,ho,HG,hg,HI,hi,spec2);
            
            bool newenv;
            if (!flag) { /*ignore this specificaiton */
                continue;
            } else { /*use this specificaiton if it has not been seen before*/
                if (ENV_HIST.size()==0) {
                    newenv = true;
                } else { /* if this is not the fist environment, then assume that the environment is seen before */
                    newenv = false;
                    for (int i=0; i<ENV_HIST.size(); i++) {
                        /* check similarity in the obstacles in the finest layer */
                        if (abs_curr->Os_[numAbs-1]->symbolicSet_!=ENV_HIST[i][0]->symbolicSet_) {
                            newenv = true;
                        } else {
                            /* check similarity in the goals in the finest layer */
                            if (abs_curr->Gs_[numAbs-1]->symbolicSet_!=ENV_HIST[i][1]->symbolicSet_) {
                                newenv = true;
                            } else {
                                /* check similarity in the initial states in the finest layer */
                                if (abs_curr->X0s_[numAbs-1]->symbolicSet_!=ENV_HIST[i][2]->symbolicSet_) {
                                    newenv = true;
                                }
                            }
                        }
                        
                        if (!newenv) /* similarity found */
                            break;
                        else if (i<ENV_HIST.size()-1) /* if this is not the last iteration of this for loop then reset newenv flag and continue with the next one */
                            newenv = false;
                    }
                }
                
                /* if the envrionment is seen before, continue with the next one */
                if (!newenv) {
                    if (verbose>0)
                        cout << "Environment seen before. Coninuing with the next one.\n";
                    continue;
                } else {
#pragma omp critical (unique_env_counter_inc)
                    unique_env_count++;
                }
            }
                      
            
            /* store the environment info in the list */
            SymbolicSet* O = new SymbolicSet(*abs_curr->Os_[numAbs-1]);
            SymbolicSet* G = new SymbolicSet(*abs_curr->Gs_[numAbs-1]);
            SymbolicSet* X0 = new SymbolicSet(*abs_curr->X0s_[numAbs-1]);
            std::array<SymbolicSet*,3> arr = { O, G, X0};
#pragma omp critical (push_new_env)
            ENV_HIST.push_back(arr);
            
            /* perform a lazy reach-avoid abstraction refinement step */
            abs_curr->onTheFlyReach(p, sys_post, radius_post, x, u);
            // debug
            abs_curr->saveFinalResult();
            //debug end
            /* print result of the abstraction refinement */
            if (abs_curr->isInitWinning()) {
                /* increment the success rate counter */
                act_success_count[1]++;
                if (useColors)
                /* print in green (the code 32) */
                    cout << "\033[32mEnvironment #" << e << "\033[0m\n";
                if (verbose>0)
                    cout << "The environment is now winning.\n";
            } else {
                if (useColors)
                    cout << "Environment #" << e << "\n";
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
            /* clear the environments and the controllers */
            abs_curr->clear_env();
            abs_curr->clear_control();
            /* update the abstraction transition relations */
#pragma omp critical (update_transitions)
            {
            for (int i=0; i<*abs_curr->system_->numAbs_; i++) {
                abs->Ts_[i]->symbolicSet_|=abs_curr->Ts_[i]->symbolicSet_;
                abs->TTs_[i]->symbolicSet_|=abs_curr->TTs_[i]->symbolicSet_;
                abs->computedDs_[i]->symbolicSet_|=abs_curr->computedDs_[i]->symbolicSet_;
            }
            } /* end of critical section*/
        } /* End of for loop over the set of environments */
        } /* end of omp parallel block */
//        delete[] abs_copies;
        if (unique_env_count==0) {
            cout << "\nSPEC value "<< spec <<" is too high. The multi-layered abstraction is too coarse. Consider recomputation with 1 more layer at the bottom.";
            return (-1);
        }
        if ((double)act_success_count[1]/unique_env_count > reqd_success_rate[1]) {
            cout << "\nTarget precision reached.";
//            break;
            /* the computation is finished, provided the SPEC value doesn't increase */
            reqd_success_rate_reached = true;
        }
        /*************************************************************/
        /* **** Computation of SPEC with abstraction refinement **** */
        /*************************************************************/
        cout << "\033[1;4;34mStarting computation of SPEC.\n\n\033[0m";
//        spec_with_no_refinement=true;
        /* iterate over all the environments */
        for (int e=0; e<NN; e++) {
            /* spawn environment with 0 distance from the given specification */
            if (!useColors)
                cout << "Environment #" << e << "\n";
            spawnO(HO,ho,verbose);
            spawnG(HG,hg,verbose);
            spawnI(HI,hi,verbose);
            
            /* *** synthesize controller on the available abstraction only *** */
            /* initialize the environment in the abstraction */
            bool flag = abs->initializeSpec(HO,ho,HG,hg,HI,hi);
            if (!flag) { /* ignore this environment */
                abs->clear_env();
                continue;
            }
            if (verbose>0)
                cout << "Abstraction initialized with the specification.\n";
            /* do a simple multi-layered reachability (without further refinement) */
            abs->plainReach(p);
            //debug
            abs->saveFinalResult();
            //debug end
            if (verbose>0)
                cout << "Controller synthesis done.\n";
            /* if the abstract game is not winning, the distance remains 0 */
            distance = 0;
            /* check if there is a controller for the abstraction */
            if (abs->isInitWinning()) { /* need to check if the controller works for the system as well */
                if (verbose>0)
                    cout << "There is a controller.\n";
                /* simulate the concrete system using the synthesized controller */
                // debug
//                abs->writeVecToFile(abs_log.trajectory,"Figures/abs_traj.txt");
                // debug end
                /* pick a random initial point within the provided initial set which is outside the obstacles and the exclusion regions in the finest layer*/
                std::vector<double> init;
                while (1) {
                    generateInitial(init,hi);
                    if (abs->X0s_[*abs->system_->numAbs_-1]->isElement(init)) {
                        break;
                    } else {
                        init.clear();
                    }
                }
                
                
                sys_log.trajectory.push_back(init);
                simulateSystem(abs,ho,hg,hi,sys_post,x,u,unsafeAt,sys_log);
                //debug
                abs->saveFinalResult();
                abs->writeVecToFile(sys_log.trajectory,"Figures/sys_traj.txt","clean");
//                abs->saveFinalResult();
                //debug end
//                printArray(sys_traj, )
                /* check if the system trajectory was unsafe */
                if (unsafeAt.size()!=0) {
                    if (useColors) {
                        /* print in red (the code 31) */
                        cout << "\033[31mEnvironment #" << e << "\033[0m\n";
                    }
                    /* recursively refine the abstraction */
                    distance = refine_recurse(abs,ho,hg,hi,sys_post,radius_post,x,u,unsafeAt,sys_log,abs_log,spec,explRadius,explHorizon,p);
                    
                    /* if spec is greater than spec_max, discontinue the process */
                    if (distance>spec_max) {
                        return -1;
                    } else if (distance>spec) {
                        spec = distance;
                    }
                    
                } else { /* if (unsafeAt.size()==0) */
                    if (useColors) {
                        /* print in green (the code 32) */
                        cout << "\033[32mEnvironment #" << e << "\033[0m\n";
                    }
                    if (verbose>0)
                        cout << "The controller worked for the system as well.\n";
                } /* End of if-else-block (unsafeAt.size()!=0) */
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
            abs->clear_env();
            abs->clear_control();
            /* clear the sets used for trajectory simulation */
            unsafeAt.clear();
            sys_log.abstraction_used.clear();
            sys_log.trajectory.clear();
            sys_log.strategy.clear();
            abs_log.abstraction_used.clear();
            abs_log.trajectory.clear();
            abs_log.strategy.clear();
        } /* End of for loop over the set of environments */
        cout << "\n\nSPEC = " << spec << "\n\n";
        
        /* if the last computed spec value worked well enough in terms of abstract game solving, and the spec value didn't increase from the last computed spec value, and finally if the abstraction was not refined further during the computation of spec, then we are done */
        if ((reqd_success_rate_reached) &&
            (spec <= spec_old) &&
            (act_success_count[0]/NN < reqd_success_rate[0])) {
//            (spec_with_no_refinement)) {
            cout << "Abstract computation finished.\nTermination condition used: reqd. abstract synthesis success rate reached + the SPEC value didn't increase afterwards + no refinement was performed during the last SPEC loop.\n";
            /* write outputs to files */
            checkMakeDir("T");
            saveVec(abs->Ts_, "T/T");
            cout << "\nWrote transitions to files. Exiting process.\n\n";
            break;
        }
        
        //debug
        /* write outputs to files */
//        checkMakeDir("T");
//        saveVec(abs->Ts_, "T/T");
        //debug end
        

        if (abs->Ts_[numAbs-1]->symbolicSet_==abs->ddmgr_->bddOne()) {
            cout << "\nNo more refinement possible.";
            /* write outputs to files */
            checkMakeDir("T");
            saveVec(abs->Ts_, "T/T");
            cout << "\nWrote transitions to files. Exiting process.\n\n";
            break;
        }
        /* clear the controller related vaiables */
        abs->clear_env();
        abs->clear_control();
//        BlackBoxReach* abs_ref = new BlackBoxReach(*abs);
        /* store the current spec value variable for later comparison */
        spec_old = spec;
//        spec = 0;
        
        iter++;
//        for (int l=0; l<dimX; l++)
//        act_success_count = 0;
    } /* End of the big while loop over spec computation and game solving */
    return spec;
}

/*! Recursively refine the abstraction for a given specification */
template<std::size_t SIZE_o, std::size_t SIZE_g, std::size_t SIZE_i, class sys_type, class rad_type, class X_type, class U_type, class L1, class L2>
double refine_recurse(BlackBoxReach* abs,
                      const std::vector<std::array<double,SIZE_o>> ho,
                      const std::vector<std::array<double,SIZE_g>> hg,
                      const std::vector<std::array<double,SIZE_i>> hi,
                      const sys_type sys_post, const rad_type radius_post,
                      const X_type x,
                      const U_type u,
                      std::vector<double>& unsafeAt,
                      L1& sys_log, L2& abs_log,
                      const double spec,
                      const X_type explRadius, const double explHorizon, const int p) {
    /* retrieve certain variables from the BlackBoxReach object */
    int verbose = abs->verbose_;
    int dimX = *abs->system_->dimX_;
    if (verbose>0) {
        cout << "System trajectory went to unsafe part at (";
        for (int i=0; i<dimX; i++) {
            cout << unsafeAt[i] << ",";
        }
        cout << "\b).\n";
    }
    /* initialize the return variable */
    double distance;
    /* save the transition BDDs of abstraction */
    std::vector<SymbolicSet*> T_old = abs->Ts_;
    /* add the problematic system trajectory to the abstraction */
    abs->addTrajectory(sys_log, explHorizon);
    cout << "\tProblematic trajectory now part of abstraction.\n";
    bool new_addition = false;
    for (int i=0; i<T_old.size(); i++) {
        if (abs->Ts_[i]->symbolicSet_!=T_old[i]->symbolicSet_) {
            new_addition = true;
            break;
        }
    }
    /* explore around the unsafe state */
//    abs->clear();
    bool new_exploration = abs->exploreAroundPoint(unsafeAt, explRadius, u, sys_post, radius_post);
    if (!new_exploration) { /* already explored upto the finest layer */
        cout << "\tNo more refinement possible.\n";
//        return;
    }
//        /* not sure if needed */
//        act_success_count[0]++;
        /* do a simple multi-layered reachability (without further refinement) */
//        abs->initializeSpec(HO,ho,HG,hg,HI,hi);
    abs->clear_control();
    abs->plainReach(p);
    if (!abs->isInitWinning()) { /* no controller exists for the abstraction */
        cout << "\tUpdated distance = 0\n";
        distance = 0;
        return distance;
    } else {
        /* there is a controller for the abstraction: system needs to be checked */
        unsafeAt.clear();
        /* simulate the system using the synthesized controller */
        /* use the old (problematic) initial state */
        std::vector<double> old_init = sys_log.trajectory[0];
        sys_log.abstraction_used.clear();
        sys_log.trajectory.clear();
        sys_log.strategy.clear();
        sys_log.trajectory.push_back(old_init);
        simulateSystem(abs,ho,hg,hi,sys_post,x,u,unsafeAt,sys_log);
        //debug
        abs->writeVecToFile(sys_log.trajectory,"Figures/sys_traj.txt","clean");
        //debug end
        if (unsafeAt.size()==0) {
            if (verbose>0)
                cout << "\tThe controller worked for the system as well.\n";
            cout << "\tUpdated distance = 0\n";
            distance = 0;
            return distance;
        }
    }
    /* if the controller failed on the system, or no more exploration was possible, recompute distance and continue refining the abstraction around the point of failure, which is stored in the variable unsafeAt */
    abs_log.abstraction_used.clear();
    abs_log.trajectory.clear();
    abs_log.strategy.clear();
    abs_log.trajectory.push_back(sys_log.trajectory[0]);
    std::vector<std::vector<double>> hovec;
    vecArr2vecVec(ho,hovec);
    abs->simulateAbs(abs_log,hovec,distance);
    cout << "\tUpdated distance = " << distance << "\n";
    //debug
    abs->writeVecToFile(abs_log.trajectory,"Figures/abs_traj.txt","clean");
    //debug end
    /* if the distance after refinement became smaller than the current spec value, no need to refine further */
    if (!new_exploration && !new_addition) {
        return distance;
    } else {
        return refine_recurse(abs,ho,hg,hi,sys_post,radius_post,x,u,unsafeAt,sys_log,abs_log,spec,explRadius,explHorizon,p);
    }
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
    int no_cont_found_count = 0;
    //debug
//    std::vector<std::vector<double>> abs_traj;
    closed_loop_log abs_log;
    double distance;
    //debug end
    for (int e=0; e<num_tests; e++) {
        bool validEnv = false;
        while (!validEnv) {
            /* clear the controller and environment info */
            abs->clear_env();
            spawnO(HO,ho,verbose);
            spawnG(HG,hg,verbose);
            spawnI(HI,hi,verbose);
            validEnv = abs->initializeSpec(HO,ho,HG,hg,HI,hi,spec_final);
        }
        if (verbose>0)
            cout << "Abstraction initialized with the specification.\n";
        /* do a simple multi-layered reachability (without further refinement) */
        abs->plainReach(p);
        if (verbose>0)
            cout << "Controller synthesis done.\n";
        if (abs->isInitWinning()) {
            if (verbose>0)
                cout << "There is a winning controller.\n";
            clog << "\nFinal number of controllers: " << abs->finalCs_.size() << '\n';
            /* now try the controller on the system itself */
            closed_loop_log sys_log, abs_log;
//            std::vector<std::vector<double>> sys_traj;
            /* pick a random initial point within the provided initial set which is outside the obstacles and the exclusion regions in the finest layer*/
            std::vector<double> init;
            while (1) {
                double toss1, toss2;
                toss1 = 0.01*(rand() % (int)(100*(hi[0][1]+hi[0][0])))-hi[0][0];
                toss2 = 0.01*(rand() % (int)(100*(hi[0][3]+hi[0][2])))-hi[0][2];
                init.push_back(toss1);
                init.push_back(toss2);
                if (abs->X0s_[*abs->system_->numAbs_-1]->isElement(init)) {
                    break;
                }
                init.clear();
            }
            
            sys_log.trajectory.push_back(init);
            std::vector<double> unsafeAt;
            simulateSystem(abs,ho,hg,hi,sys_post,x,u,unsafeAt,sys_log);
            // debug
            abs->saveFinalResult(); abs->writeVecToFile(sys_log.trajectory,"Figures/sys_traj.txt","clean");
            abs_log.trajectory.clear();
            abs_log.strategy.clear();
            abs_log.trajectory.push_back(sys_log.trajectory[0]);
            std::vector<std::vector<double>> hovec;
            vecArr2vecVec(ho,hovec);
            abs->simulateAbs(abs_log,hovec,distance);
            abs->writeVecToFile(abs_log.trajectory,"Figures/abs_traj.txt","clean");
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
            no_cont_found_count++;
        }
        // debug
        abs->saveFinalResult();
        //end
        /* clear the specification sets */
        ho.clear();
        hg.clear();
        hi.clear();
        HO.clear();
        HG.clear();
        HI.clear();
        /* clear the controller and the environment */
        abs->clear_env();
        abs->clear_control();
        
        
        //        /* save synthesis results */
        //        abs->saveFinalResult();
    }
    cout << "\nTest result:\n";
    if (num_tests!=no_cont_found_count) {
        cout << "Success rate = " << success_count*100/(num_tests-no_cont_found_count) << "% (" << success_count << " times).\n";
        if (no_cont_found_count!=0) {
            cout << "This excludes the " << no_cont_found_count*100/num_tests << "% (" << no_cont_found_count << ") cases where no controller could be synthesized.";
        }
    }
    else
        cout << "All tests failed.\n";
}

#endif
