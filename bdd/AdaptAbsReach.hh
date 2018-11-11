#ifndef ADAPTABSREACH_HH_
#define ADAPTABSREACH_HH_

#include <cstdio>
#include <vector>

#include "SymbolicModelGrowthBound.hh"
#include "System.hh"
#include "Helper.hh"

using std::clog;
using std::freopen;
using std::string;
using std::vector;
using namespace helper;

namespace scots {

enum ReachResult {CONVERGEDVALID, CONVERGEDINVALID, NOTCONVERGED};

class AdaptAbsReach {
public:
    Cudd* ddmgr_; /*!< A single manager object common to all BDDs used in the program. */
    System* system_; /*!< Contains abstraction parameters. */
    vector<double*> etaXs_; /*!< *system_->numAbs_ x *system_->dimX_ matrix of state space grid spacings. */
    vector<double*> tau_; /*!< *system_->numAbs_ x 1 matrix of time steps. */
    vector<SymbolicSet*> Xs_; /*!< The *system_->numAbs_ "pre" state space abstractions, coarsest (0) to finest. */
    vector<SymbolicSet*> Os_; /*!< Instance of *Xs_[i] containing unsafe (obstacle) states. */
    vector<SymbolicSet*> Zs_; /*!< Instance of *Xs_[i] containing winning states. */
    vector<SymbolicSet*> X2s_; /*!< The *system_->numAbs_ "post" state space abstractions, coarsest (0) to finest. */
    SymbolicSet* U_; /*!< The single input space abstraction. */
    vector<SymbolicSet*> Cs_; /*!< Controller \subseteq *Xs_[i] x *U_. */
    int numBDDVars_; /*!< Total number of BDD variables used by ddmgr_. */
    vector<BDD*> cubesX_; /*!< *cubesX_[i] is a BDD with a single minterm of all 1s over the domain of *Xs_[i]. */
    vector<BDD*> cubesX2_; /*!< *cubesX2_[i] is a BDD with a single minterm of all 1s over the domain of *X2s_[i]. */
    BDD* cubeU_; /*!< A BDD with a single minterm of all 1s over the domain of the input SymbolicSet. */
    vector<BDD*> notXUvars_; /*!< notXUvars_[i] is a BDD over all ddmgr_ variables whose 1-term is DC for *Xs_[i] and *U_ and 1 otherwise. */
    vector<BDD*> notXvars_; /*!< notXvars_[i] is a BDD over all ddmgr_ variables whose 1-term is DC for *Xs_[i] and 1 otherwise. */
    vector<BDD*> cubesCoarser_; /*!< Helper BDD for projection to coarser layer. */
    vector<SymbolicSet*> Ts_; /*!< Ts_[i] stores the transition relation of abstraction i. */
    vector<SymbolicSet*> TTs_; /*!< TTs_[i] is Ts_[i] with X2s_[i] existentially abtracted. */
    vector<int*> permutesXtoX2_; /*!< To transform a BDD over X variables to an equivalent one over X2 variables. */
    vector<int*> permutesX2toX_; /*!< To transform a BDD over X2 variables to an equivalent one over X variables. */
    vector<int*> permutesCoarser_; /*!< To project a BDD to the next coarser layer. */
    vector<int*> permutesFiner_; /*!< To project a BDD to the next finer layer. */
    vector<OdeSolver*> solvers_; /*!< ODE solvers (Runge-Katta approximation) for each abstraction time step. */

    vector<SymbolicSet*> Gs_; /*!< Instance of *Xs_[i] containing goal states. */
    vector<SymbolicSet*> validZs_; /*!< Contains winning states that act as savepoints. */
    vector<SymbolicSet*> validCs_; /*!< Controllers that act as savepoints. */
    vector<SymbolicSet*> finalCs_; /*!< Sequence of controllers that satisfy the specification. */
    vector<SymbolicSet*> finalZs_; /*!< Sequence of domains of finalCs_. */
    vector<int> finalAbs_; /*!< Sequence of abstraction layer indices of finalCs_. */
    vector<SymbolicSet*> Ds_; /*!< Instance of *Xs_[i] containing possible winning states (outer projections). */
    vector<SymbolicSet*> innerDs_; /*!< Instance of *Xs_[i] containing possible winning states (inner projections). */
    vector<SymbolicSet*> computedDs_; /*!< Instance of *Xs_[i] containing states for which transitions have already been constructed. */

    vector<SymbolicSet*> uTs_; /*!< Transition relations for exploring, i.e. all with coarsest space gridding but varying time sampling. */

    int m_; /*!< Number of iterations for synthesis in a layer before attempting to go coarser. */
    int p_; /*!< Number of iterations for exploration before attempting synthesis. */

    double abTime_; /*!< Abstraction time. */
    double synTime_; /*!< Synthesis time. */
    int verbose_; /*!< 0(default)=do not print intermediate result, 1=print intermediate result on std I/O. */

    /*!	Constructor for an AdaptAbsReach object.
     *  \param[in]  logFile     Filename of program log.
     */
    AdaptAbsReach(const char* logFile, int verbose=0) {
        freopen(logFile, "w", stderr);
        clog << logFile << '\n';

        abTime_ = 0;
        synTime_ = 0;
        verbose_ = verbose;
    }
    /*! Destructor for an AdaptAbsReach object. */
    ~AdaptAbsReach() {
        clog << "Abstraction construction time: " << abTime_ << " seconds.\n";
        clog << "Controller synthesis time: " << synTime_ << " seconds.\n";
        clog << "Abs + syn time: " << abTime_ + synTime_ << " seconds.\n";
        deleteVecArray(etaXs_);
        deleteVec(tau_);
        deleteVec(Xs_);
        deleteVec(Os_);
        deleteVec(Zs_);
        deleteVec(X2s_);
        delete U_;
        deleteVec(Cs_);
        deleteVec(cubesX_);
        deleteVec(cubesX2_);
        delete cubeU_;
        deleteVec(notXUvars_);
        deleteVec(notXvars_);
        deleteVec(cubesCoarser_);
        deleteVec(Ts_);
        deleteVec(TTs_);
        deleteVecArray(permutesXtoX2_);
        deleteVecArray(permutesX2toX_);
        deleteVecArray(permutesCoarser_);
        deleteVecArray(permutesFiner_);
        deleteVec(solvers_);
        deleteVec(Gs_);
        deleteVec(validZs_);
        deleteVec(validCs_);
        deleteVec(finalCs_);
        deleteVec(finalZs_);
        deleteVec(Ds_);
        deleteVec(innerDs_);
        deleteVec(computedDs_);
        deleteVec(uTs_);
        fclose(stderr);
        delete ddmgr_;
    }

    /*! Non-lazy reachability wrapper function.
    *  \param[in]  p           Parameter controlling amount of finer layer synthesis attempts.
    *  \param[in]  readAbs     Flag set to TRUE when abstractions are to be read from file.
    *  \param[in]  sysNext     System ODE.
    *  \param[in]  radNext     Growth bound ODE.
    *  \param[in]  x           Dummy state point.
    *  \param[in]  u           Dummy input point.
    */
    template<class sys_type, class rad_type, class X_type, class U_type>
    void eagerReach(int p, bool readAbs, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        m_ = p; // max. iterations for consecutive reachability for non-coarsest layers
        clog << "m: " << m_ << '\n';

        // pre-compute all abstractions
        TicToc timer;
        timer.tic();
        if (!readAbs) {
            for (int ab = 0; ab < *system_->numAbs_; ab++) {
                Ds_[ab]->addGridPoints();
                computeAbstraction(ab, sysNext, radNext, x, u);
            }
            clog << "Ts_ computed by sampling system behavior.\n";
        }
        else {
            loadTs();
        }

        abTime_ += timer.toc();

        if (!readAbs) {
            checkMakeDir("T");
            saveVec(Ts_, "T/T");
            clog << "Wrote Ts_ to file.\n";
        }

        // synthesize controller
        int ab = 0;
        eagerReachRecurse(ab);

        clog << "\nFinal number of controllers: " << finalCs_.size() << '\n';

        checkMakeDir("C");
        saveVec(finalCs_, "C/C");
        checkMakeDir("Z");
        saveVec(finalZs_, "Z/Z");
        checkMakeDir("T");
        saveVec(Ts_, "T/T");
        clog << "Wrote Ts_ to file.\n";
        return;
    }

    /*! Recursive synthesis of non-lazy controller.
     *  \param[in]  ab          Current abstraction.
     */
     void eagerReachRecurse(int ab) {
         clog << '\n';
         clog << "current abstraction: " << ab << '\n';
         clog << "controllers: " << finalCs_.size() << '\n';
         ReachResult result;
         if (ab == 0) {
             result = reach(ab);
         }
         else {
             result = reach(ab, m_);
         }
         if (result == CONVERGEDVALID) {
             clog << "result: converged valid\n";
         }
         else if (result == CONVERGEDINVALID) {
             clog << "result: converged invalid\n";
         }
         else {
             clog << "result: not converged\n";
         }

         if (result != CONVERGEDINVALID) {
             saveCZ(ab);
             validZs_[ab]->symbolicSet_ = Zs_[ab]->symbolicSet_;
             validCs_[ab]->symbolicSet_ = Cs_[ab]->symbolicSet_;
             clog << "saved as snapshot, saved to valids\n";
         }
         else { // result is CONVERGEDINVALID
             Zs_[ab]->symbolicSet_ = validZs_[ab]->symbolicSet_; // reset this layer's progress
             Cs_[ab]->symbolicSet_ = validCs_[ab]->symbolicSet_;
             clog << "reset to valids\n";
         }

         if (result != NOTCONVERGED) {
             if (ab == *system_->numAbs_ - 1) {
                 return;
             }
             else { // go finer
                 clog << "Going finer\n";
                 int nextAb = ab + 1;
                 finer(Zs_[ab], Zs_[nextAb], ab);
                 Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
                 validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;
                 eagerReachRecurse(nextAb);
                 return;
             }
         }
         else { // not converged, go coarser
             clog << "Going coarser\n";
             int nextAb = ab - 1;
             coarserInner(Zs_[nextAb], Zs_[ab], nextAb);
             Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
             validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;
             eagerReachRecurse(nextAb);
             return;
         }
     }

    /*! Lazy reachability wrapper function.
     *  \param[in]  p           Parameter controlling amount of exploration to do between synthesis attempts.
     *  \param[in]  sysNext     System ODE.
     *  \param[in]  radNext     Growth bound ODE.
     *  \param[in]  x           Dummy state point.
     *  \param[in]  u           Dummy input point.
     */
    template<class sys_type, class rad_type, class X_type, class U_type>
    void onTheFlyReach(int p, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        m_ = p; // max. iterations for consecutive reachability for non-coarsest layers
        p_ = p; // coarsest layer uncontrollable-pred. parameter
        clog << "m: " << m_ << '\n';
        clog << "p: " << p_ << '\n';

        // start by synthesizing all uTs_
        Ds_[0]->addGridPoints();

        TicToc timer;
        timer.tic();
        computeExplorationAbstractions(sysNext, radNext, x, u);
        abTime_ += timer.toc();
        cout << "\nComputation time for exploration abstractions: " << timer.toc() <<"\n";

        // Ts_[0] is uTs_[0] with obstacles removed
        Ts_[0]->symbolicSet_ = uTs_[0]->symbolicSet_ & !Os_[0]->symbolicSet_;
        TTs_[0]->symbolicSet_ = Ts_[0]->symbolicSet_.ExistAbstract(*notXUvars_[0]);

        // begin on-the-fly reachability synthesis
        int ab = 0;
        onTheFlyReachRecurse(ab, sysNext, radNext, x, u);

        clog << "\nFinal number of controllers: " << finalCs_.size() << '\n';

        checkMakeDir("C");
        saveVec(finalCs_, "C/C");
        checkMakeDir("Z");
        saveVec(finalZs_, "Z/Z");
        checkMakeDir("T");
        saveVec(Ts_, "T/T");
        clog << "Wrote Ts_ to file.\n";

        if (verbose_) {
          for (size_t i = 0; i < finalCs_.size(); i++) {
            std::cout << "Controller " << i+1 << ":"<< '\n';
            finalCs_[i]->printInfo();
            std::cout << "Target " << i+1 << '\n';
            finalZs_[i]->printInfo();
          }
        }
        return;
    }

    /*! Main lazy reachability synthesis procedure.
     *  \param[in]  ab          Layer 0-index.
     *  \param[in]  sysNext     System ODE.
     *  \param[in]  radNext     Growth bound ODE.
     *  \param[in]  x           Dummy state point.
     *  \param[in]  u           Dummy input point.
     */
    template<class sys_type, class rad_type, class X_type, class U_type>
    void onTheFlyReachRecurse(int ab, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        clog << '\n';
        clog << "current abstraction: " << ab << '\n';
        clog << "controllers: " << finalCs_.size() << '\n';

        if (verbose_) {
            cout << '\n';
            cout << "current abstraction: " << ab << '\n';
            cout << "controllers: " << finalCs_.size() << '\n';
        }

        ReachResult result;
        if (ab == 0) {
            TicToc timer;
            timer.tic();
            result = reach(ab);
            synTime_ += timer.toc();
        }
        else {
            TicToc timer;
            timer.tic();
            result = reach(ab, m_);
            synTime_ += timer.toc();
        }
        if (result == CONVERGEDVALID) {
            clog << "result: converged valid\n";
            if (verbose_)
                cout << "result: converged valid\n";
        }
        else if (result == CONVERGEDINVALID) {
            clog << "result: converged invalid\n";
            if (verbose_)
                cout << "result: converged invalid\n";
        }
        else {
            clog << "result: not converged\n";
            if (verbose_)
                cout << "result: not converged\n";
        }

        if (result != CONVERGEDINVALID) {
            saveCZ(ab);
            validZs_[ab]->symbolicSet_ = Zs_[ab]->symbolicSet_;
            validCs_[ab]->symbolicSet_ = Cs_[ab]->symbolicSet_;
            clog << "saved as snapshot, saved to valids\n";
            if (verbose_)
                cout << "saved as snapshot, saved to valids\n";
        }
        else { // result is CONVERGEDINVALID
            Zs_[ab]->symbolicSet_ = validZs_[ab]->symbolicSet_; // reset this layer's progress
            Cs_[ab]->symbolicSet_ = validCs_[ab]->symbolicSet_;
            clog << "reset to valids\n";
            if (verbose_)
                cout << "reset to valids\n";
        }

        if (result != NOTCONVERGED) { // ab = 0 always converges
            if (ab == *system_->numAbs_ - 1) {
                return;
            }
            else { // go finer
                clog << "Going finer\n";
                if (verbose_)
                    cout << "Going finer\n";
                int nextAb = ab + 1;
                TicToc timer;
                timer.tic();
                explore(ab, nextAb, sysNext, radNext, x, u);
                abTime_ += timer.toc();
                finer(Zs_[ab], Zs_[nextAb], ab);
                Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
                validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;
                onTheFlyReachRecurse(nextAb, sysNext, radNext, x, u);
                return;
            }
        }
        else { // not converged, go coarser
            clog << "Going coarser\n";
            if (verbose_)
                cout << "Going coarser\n";
            int nextAb = ab - 1;
            if (nextAb != 0) { // pointless to do for coarsest layer
                TicToc timer;
                timer.tic();
                explore(ab, nextAb, sysNext, radNext, x, u);
                abTime_ += timer.toc();
            }
            coarserInner(Zs_[nextAb], Zs_[ab], nextAb);
            Zs_[nextAb]->symbolicSet_ |= validZs_[nextAb]->symbolicSet_;
            validZs_[nextAb]->symbolicSet_ = Zs_[nextAb]->symbolicSet_;
            onTheFlyReachRecurse(nextAb, sysNext, radNext, x, u);
            return;
        }
    }

    /*! Fully computes the exploration transition relations (which use the coarsest state gridding).
     *  \param[in]  sysNext     System ODE.
     *  \param[in]  radNext     Growth bound ODE.
     *  \param[in]  x           Dummy state point.
     *  \param[in]  u           Dummy input point.
     */
    template<class sys_type, class rad_type, class X_type, class U_type>
    void computeExplorationAbstractions(sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        if (verbose_==1)
          cout << "Computing exploration abstractions.\n";
        SymbolicModelGrowthBound<X_type, U_type> abstraction(Ds_[0], U_, X2s_[0]); // coarsest state gridding
        for (int ab = 0; ab < *system_->numAbs_; ab++) {
            // abstraction.computeTransitionRelation(sysNext, radNext, *solvers_[ab], verbose_); // use solver with time step corresponding to that of each layer
            uTs_[ab]->symbolicSet_ = abstraction.transitionRelation_; // add to transition relation
        }
        abstraction.computeExplorationTransitionRelation(sysNext, radNext, *system_->numAbs_, solvers_, uTs_, verbose_);
        x = x; // gets rid of warning message regarding lack of use
        u = u;
    }

    /*! Does lazy exploration based on the current status of controller synthesis.
     *  \param[in]  curAb       The current layer 0-index.
     *  \param[in]  nextAb      The next layer 0-index.
     *  \param[in]  sysNext     System ODE.
     *  \param[in]  radNext     Growth bound ODE.
     *  \param[in]  x           Dummy state point.
     *  \param[in]  u           Dummy input point.
     */
    template<class sys_type, class rad_type, class X_type, class U_type>
    void explore(int curAb, int nextAb, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        Ds_[curAb]->symbolicSet_ = validZs_[curAb]->symbolicSet_ | Gs_[curAb]->symbolicSet_; // target for uncontrollable reach
        innerDs_[curAb]->symbolicSet_ = validZs_[curAb]->symbolicSet_ | Gs_[curAb]->symbolicSet_;
        for (int c = curAb - 1; c >= 0; c--) {
            coarserOuter(Ds_[c], Ds_[c+1], c);
            coarserInner(innerDs_[c], innerDs_[c+1], c);
        }
        Ds_[0]->symbolicSet_ = cooperativeReach(Ds_[0]->symbolicSet_, curAb, p_); // do cooperative reach
        Ds_[0]->symbolicSet_ &= (!innerDs_[0]->symbolicSet_); // states to do abstraction for
        for (int c = 0; c < nextAb; c++) {
            finer(Ds_[c], Ds_[c+1], c);
        }
        computeAbstraction(nextAb, sysNext, radNext, x, u);
    }

    /*! Calculates the transitions stemming from a certain subset of the state space D and adds them to the transition relation.
     *  \param[in]  ab      The layer 0-index.
     *  \param[in]  sysNext     System ODE.
     *  \param[in]  radNext     Growth bound ODE.
     *  \param[in]  x           Dummy state point.
     *  \param[in]  u           Dummy input point.
     */
    template<class sys_type, class rad_type, class X_type, class U_type>
    void computeAbstraction(int ab, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
        BDD D = Ds_[ab]->symbolicSet_ & (!computedDs_[ab]->symbolicSet_); // don't repeat sampling
        D &= !Os_[ab]->symbolicSet_;
        computedDs_[ab]->symbolicSet_ |= Ds_[ab]->symbolicSet_; // update computed part of transition relation
        Ds_[ab]->symbolicSet_ = D;

        SymbolicModelGrowthBound<X_type, U_type> abstraction(Ds_[ab], U_, X2s_[ab]); // abstraction computation will only iterate over the elements in the domain of the state space SymbolicSet (first argument)
        abstraction.computeTransitionRelation(sysNext, radNext, *solvers_[ab], verbose_); // was hard-coded to 0, source of "tunneling" bug
        if (abstraction.transitionRelation_ != ddmgr_->bddZero()) { // no point adding/displaying if nothing was added
            Ts_[ab]->symbolicSet_ |= abstraction.transitionRelation_; // add to transition relation
//            BDD O2 = Os_[ab]->symbolicSet_.Permute(permutesXtoX2_[ab]); // causes unsound behavior, leaving here as warning
//            Ts_[ab]->symbolicSet_ &= !O2;
//            Ts_[ab]->printInfo(1);
            TTs_[ab]->symbolicSet_ = Ts_[ab]->symbolicSet_.ExistAbstract(*notXUvars_[ab]);
        }
        x = x; // gets rid of warning message regarding lack of use
        u = u;
    }

    /*! Calculates the reachability fixed-point for m iterations.
     *  \param[in]  ab      The layer 0-index.
     *  \param[in]  m       The number of iterations. Default -1 corresponds to infinity.
     *  \returns    Integer representing status of convergence.
     */
    ReachResult reach(int ab, int m = -1) {
        int i = 1;
        clog << "abstraction: " << ab << '\n';
        if (verbose_)
            cout << "abstraction: " << ab << '\n';
        while (1) {
            clog << "iteration: " << i << '\n';
            if (verbose_)
                cout << "iteration: " << i << '\n';
            BDD controllablePreZ = controllablePre(Zs_[ab]->symbolicSet_, ab);
            BDD C = controllablePreZ | Gs_[ab]->symbolicSet_;
            BDD N = C & (!(Cs_[ab]->symbolicSet_.ExistAbstract(*cubeU_)));
            Cs_[ab]->symbolicSet_ |= N;
            Zs_[ab]->symbolicSet_ = C.ExistAbstract(*notXvars_[ab]) | validZs_[ab]->symbolicSet_; // the disjunction part is new

            if (N == ddmgr_->bddZero() && i != 1) {
                if (i >= 2) {
                    return CONVERGEDVALID;
                }
                else{
                    return CONVERGEDINVALID;
                }
            }
            if (i == m) {
                return NOTCONVERGED;
            }
            i += 1;
        }
    }

    /*! Calculates the cooperative reachability fixed-point for p iterations.
     *  \param[in]  Z       The given set.
     *  \param[in]  ab      The layer 0-index of Z.
     *  \param[in]  p       The number of iterations.
     *  \return     BDD containing all possible winning states at most p transitions away from Z.
     */
    BDD cooperativeReach(BDD Z, int ab, int p) {
        int i = 1;
        while (1) {
            BDD cooperativePreZ = cooperativePre(Z, ab);
            BDD N = cooperativePreZ & (!Z);
            Z = cooperativePreZ;

            if (i == p || N == ddmgr_->bddZero()) {
                return Z;
            }
            i += 1;
        }
    }

    /*! Calculates the cooperative predecessor of a given set with respect to the specified un-restricted transition relation.
     *  \param[in]  Z       The given set.
     *  \param[in]  ab      0-index of the un-restricted transition relation.
     *  \return     BDD containing {x} for which there exists a u s.t. there exists a post state of (x,u) in Z.
     */
    BDD cooperativePre(BDD Z, int ab) {
        BDD Z2 = Z.Permute(permutesXtoX2_[0]); // ab = 0 because we use the coarsest state griddings
        BDD uPreZ = uTs_[ab]->symbolicSet_.AndAbstract(Z2, *notXUvars_[0]);
        uPreZ = uPreZ.ExistAbstract(*notXvars_[0]);
        return uPreZ;
    }


    /*! Calculates the controllable predecessor of the given set with respect to the transition relation at the specified level of abstraction.
     *  \param[in]  Z           The winning set.
     *  \param[in]  ab      0-index of the current abstraction.
     *  \return     BDD containing {(x,u)} for which all post states are in Z.
     */
    BDD controllablePre(BDD Z, int ab) {
        // swap to X2s_[ab]
        BDD Z2 = Z.Permute(permutesXtoX2_[ab]);
        // posts outside of winning set
        BDD nZ2 = !Z2;
        // {(x,u)} with some posts outside of winning set
        BDD Fbdd = Ts_[ab]->symbolicSet_.AndAbstract(nZ2, *notXUvars_[ab]);
        // {(x,u)} with no posts outside of winning set
        BDD nF = !Fbdd;
        // get rid of junk
        BDD preF = TTs_[ab]->symbolicSet_.AndAbstract(nF, *notXUvars_[ab]);
        return preF;
    }

    /*! Calculates the finer projection of Zc and stores it in Zf.
     *  \param[in]  Zc      Contains set to project.
     *  \param[in]  Zf      Contains result of projection.
     *  \param[in]  c       Layer 0-index of Zc.
     */
    void finer(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
        Zf->addGridPoints();
        Zf->symbolicSet_ &= Zc->symbolicSet_.Permute(permutesFiner_[c]);
    }

    /*! Calculates the coarser inner projection of Zf and stores it in Zc.
     *  \param[in]  Zc      Contains result of projection.
     *  \param[in]  Zf      Contains set to project.
     *  \param[in]  c       Layer 0-index of Zc.
     */
    void coarserInner(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
        Zc->symbolicSet_ = (Zf->symbolicSet_.UnivAbstract(*cubesCoarser_[c])).Permute(permutesCoarser_[c]);
    }

    /*! Calculates the coarser outer projection of Zf and stores it in Zc.
     *  \param[in]  Zc      Contains result of projection.
     *  \param[in]  Zf      Contains set to project.
     *  \param[in]  c       Layer 0-index of Zc.
     */
    void coarserOuter(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
        Zc->symbolicSet_ = (Zf->symbolicSet_.ExistAbstract(*cubesCoarser_[c])).Permute(permutesCoarser_[c]);
    }

    /*! Initializes data members.
     *  \param[in]  system      Contains abstraction parameters.
     *  \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
     *  \param[in]  addG    Function pointer specifying the points that should be added to the goal set.
     */
    template<class O_type, class G_type>
    void initialize(System* system, O_type addO, G_type addG) {
        ddmgr_ = new Cudd;
        system_ = system;

        TicToc timer;
        timer.tic();

        // This order should be kept.
        initializeEtaTau();
        initializeSolvers();
        initializeSymbolicSets(addO, addG);
        initializeNumBDDVars();
        initializePermutes();
        initializeCubes();
        initializeNotVars();
        initializeProjs();
        verifySave();

        abTime_ += timer.toc();
    }

    /*! Initializes SymbolicSet data members.
     *  \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
     *  \param[in]	addG	Function pointer specifying the points that should be added to the goal set.
     */
    template<class O_type, class G_type>
    void initializeSymbolicSets(O_type addO, G_type addG) {
        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* X = new SymbolicSet(*ddmgr_, *system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[i], tau_[i][0]);
            X->addGridPoints();
            Xs_.push_back(X);
        }
        clog << "Xs_ initialized with full domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* X2 = new SymbolicSet(*Xs_[i], 1);
            X2->addGridPoints();
            X2s_.push_back(X2);
        }
        clog << "X2s_ initialized with full domain.\n";

        U_ = new SymbolicSet(*ddmgr_, *system_->dimU_, system_->lbU_, system_->ubU_, system_->etaU_, 0);
        U_->addGridPoints();
        clog << "U_ initialized with full domain.\n";

        initializeSpec(&Os_, addO);
        clog << "Os_ initialized according to specification.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* Z = new SymbolicSet(*Xs_[i]);
            Zs_.push_back(Z);

        }
        clog << "Zs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* C = new SymbolicSet(*Xs_[i], *U_);
            Cs_.push_back(C);
        }
        clog << "Cs_ initialized with empty domain.\n";

        initializeSpec(&Gs_, addG);
        clog << "Gs_ initialized according to specification.\n";

        // checking that specification is valid
        for (int i = 0; i < *system_->numAbs_; i++) {
            if ((Gs_[i]->symbolicSet_ & Os_[i]->symbolicSet_) != ddmgr_->bddZero()) {
                error("Error: G and O have nonzero intersection.\n");
            }
        }
        clog << "No obstacle problem with specification.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* validZ = new SymbolicSet(*Xs_[i]);
            validZs_.push_back(validZ);
        }
        clog << "validZs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* validC = new SymbolicSet(*Xs_[i], *U_);
            validCs_.push_back(validC);
        }
        clog << "validCs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* D = new SymbolicSet(*Xs_[i]);
            Ds_.push_back(D);
        }
        clog << "Ds_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* innerD = new SymbolicSet(*Xs_[i]);
            innerDs_.push_back(innerD);
        }
        clog << "innerDs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* computedD = new SymbolicSet(*Xs_[i]);
            computedDs_.push_back(computedD);
        }
        clog << "computedDs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* T = new SymbolicSet(*Cs_[i], *X2s_[i]);
            Ts_.push_back(T);
        }
        clog << "Ts_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* TT = new SymbolicSet(*Cs_[i]);
            TTs_.push_back(TT);
        }
        clog << "TTs_ initialized with empty domain.\n";

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* uT = new SymbolicSet(*Ts_[0]); // exploration transition relations are all with coarsest gridding
            uTs_.push_back(uT);
        }
        clog << "uTs_ initialized with empty domain.\n";
    }

    /*! Initializes a vector of SymbolicSets that is an instance of Xs_ containing points as specified by addSpec.
     *  \param[in] vec      Vector of SymbolicSets to initialize.
     *  \param[in] addSpec  Function pointer specifying the points to add to each element of vec.
     */
    template<class vec_type, class spec_type>
    void initializeSpec(vec_type* vec, spec_type addSpec) {
        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* Xinstance = new SymbolicSet(*Xs_[i]);
            addSpec(Xinstance);
            vec->push_back(Xinstance);
        }
    }

    /*! Initializes the BDDs useful for existential abstraction. */
    void initializeNotVars() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            BDD* notXUvars = new BDD;
            BDD* notXvars = new BDD;
            *notXUvars = ddmgr_->bddOne();
            for (int j = 0; j < *system_->numAbs_; j++) {
                *notXUvars &= *cubesX2_[j];
                if (i != j) {
                    *notXUvars &= *cubesX_[j];
                }
            }
            *notXvars = *notXUvars & *cubeU_;

            notXUvars_.push_back(notXUvars);
            notXvars_.push_back(notXvars);
        }

        clog << "Initialized notVars.\n";
    }

    /*! Initializes the abstractions' state space grid parameters and time sampling parameters. */
    void initializeEtaTau() {
        for (int i = 0; i < *system_->dimX_; i++) {
            if (system_->etaRatio_[i] != 1 && system_->etaRatio_[i] != 2) {
                error("Error: unsupported etaRatio, which must consist of only 1s or 2s.\n");
            }
        }

        double* etaCur = new double[*system_->dimX_];
        for (int i = 0; i < *system_->dimX_; i++) {
            etaCur[i] = system_->etaX_[i];
        }
        double* tauCur = new double;
        *tauCur = *system_->tau_;

        for (int i = 0; i < *system_->numAbs_; i++) {
            double* etaX = new double[*system_->dimX_];
            double* taui = new double;
            for (int j = 0; j < *system_->dimX_; j++) {
                etaX[j] = etaCur[j];
            }
            *taui = *tauCur;
            etaXs_.push_back(etaX);
            tau_.push_back(taui);

            for (int j = 0; j < *system_->dimX_; j++) {
                etaCur[j] /= system_->etaRatio_[j];
            }
            *tauCur /= *system_->tauRatio_;
        }

        clog << "Initialized etaXs_, tau_.\n";

        delete[] etaCur;
        delete tauCur;
    }

    /*! Initializes the Runge-Katta ODE solvers.
     */
    void initializeSolvers() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            OdeSolver* solver = new OdeSolver(*system_->dimX_, *system_->nSubInt_, *tau_[i]);
            solvers_.push_back(solver);
        }
        clog << "Initialized solvers.\n";
    }

    /*! Initializes the number of total distinct BDD variables in use. */
    void initializeNumBDDVars() {
        numBDDVars_ = 0;
        for (int i = 0; i < *system_->numAbs_; i++) {
            numBDDVars_ += Xs_[i]->nvars_;
            numBDDVars_ += X2s_[i]->nvars_;
        }
        numBDDVars_ += U_->nvars_;

        clog << "Number of BDD variables: " << numBDDVars_ << '\n';
    }

    /*! Initializes the BDDs used for efficient projection between consecutive layers of abstraction. */
    void initializeProjs() {
        int ones = 0;
        for (int i = 0; i < *system_->dimX_; i++) {
            if (system_->etaRatio_[i] == 1) {
                ones++;
            }
        }

        for (int i = 1; i < *system_->numAbs_; i++) {
            BDD* varsCoarser = new BDD[*system_->dimX_ - ones];
            int ind = 0;
            for (int j = 0; j < *system_->dimX_; j++) {
                if (system_->etaRatio_[j] == 2) {
                    varsCoarser[ind] = ddmgr_->bddVar(Xs_[i]->indBddVars_[j][0]);
                    ind++;
                }
            }
            BDD* cubeCoarser = new BDD;
            *cubeCoarser = ddmgr_->bddComputeCube(varsCoarser, NULL, *system_->dimX_ - ones);
            cubesCoarser_.push_back(cubeCoarser);
            delete[] varsCoarser;
        }

        for (int i = 1; i < *system_->numAbs_; i++) {
            int* permuteCoarser = new int[numBDDVars_];
            for (int ind = 0; ind < numBDDVars_; ind++) {
                permuteCoarser[ind] = 0;
            }

            for (int dim = 0; dim < *system_->dimX_; dim++) {
                if (system_->etaRatio_[dim] == 2) {
                    for (size_t projInd = 1; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                        permuteCoarser[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i-1]->indBddVars_[dim][projInd-1];
                    }
                }
                else if (system_->etaRatio_[dim] == 1){
                    for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                        permuteCoarser[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i-1]->indBddVars_[dim][projInd];
                    }
                }
            }
            permutesCoarser_.push_back(permuteCoarser);

            clog << "permuteCoarser " << i << " to " << i-1 << ": ";
            printArray(permuteCoarser, numBDDVars_);
        }

        for (int i = 0; i < *system_->numAbs_ - 1; i++) {
            int* permuteFiner = new int[numBDDVars_];
            for (int ind = 0; ind < numBDDVars_; ind++) {
                permuteFiner[ind] = 0;
            }

            for (int dim = 0; dim < *system_->dimX_; dim++) {
                if (system_->etaRatio_[dim] == 2) {
                    for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                        permuteFiner[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i+1]->indBddVars_[dim][projInd+1];
                    }
                }
                else if (system_->etaRatio_[dim] == 1) {
                    for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                        permuteFiner[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i+1]->indBddVars_[dim][projInd];
                    }
                }
            }
            permutesFiner_.push_back(permuteFiner);

            clog << "permuteFiner " << i << " to " << i+1 << ": ";
            printArray(permuteFiner, numBDDVars_);
        }
    }

    /*! Initializes the BDDs that are the precursors to notXvars_ and notXUvars_ */
    void initializeCubes() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            BDD* cubeX = new BDD;
            BDD* cubeX2 = new BDD;

            BDD* varsX = new BDD[Xs_[i]->nvars_];
            BDD* varsX2 = new BDD[X2s_[i]->nvars_];

            for (size_t j = 0; j < Xs_[i]->nvars_; j++) {
                varsX[j] = ddmgr_->bddVar(Xs_[i]->idBddVars_[j]);
            }

            for (size_t j = 0; j < X2s_[i]->nvars_; j++) {
                varsX2[j] = ddmgr_->bddVar(X2s_[i]->idBddVars_[j]);
            }

            *cubeX = ddmgr_->bddComputeCube(varsX, NULL, Xs_[i]->nvars_);
            *cubeX2 = ddmgr_->bddComputeCube(varsX2, NULL, X2s_[i]->nvars_);


            cubesX_.push_back(cubeX);
            cubesX2_.push_back(cubeX2);

            delete[] varsX;
            delete[] varsX2;
        }

        cubeU_ = new BDD;
        *cubeU_ = U_->getCube();

        clog << "Initialized cubes.\n";
    }

    /*! Initializes the arrays of BDD variable IDs that allow a BDD over a pre/post-X domain to be projected to the identical BDD over the corresponding post/pre-X domain.
     */
    void initializePermutes() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            int* permuteXtoX2 = new int[numBDDVars_];
            int* permuteX2toX = new int[numBDDVars_];
            for (int j = 0; j < numBDDVars_; j++) {
                permuteXtoX2[j] = 0;
                permuteX2toX[j] = 0;
            }

            for (size_t j = 0; j < Xs_[i]->nvars_; j++) {
                permuteXtoX2[Xs_[i]->idBddVars_[j]] = X2s_[i]->idBddVars_[j];
                permuteX2toX[X2s_[i]->idBddVars_[j]] = Xs_[i]->idBddVars_[j];
            }

            permutesXtoX2_.push_back(permuteXtoX2);
            permutesX2toX_.push_back(permuteX2toX);
        }

        clog << "Initialized permutes.\n";
    }

    /*! Saves and prints to log file some information related to the reachability specification. */
    void verifySave() {
        clog << "etaXs_:\n";
        for (size_t i = 0; i < etaXs_.size(); i++) {
            clog << "abstraction " << i << ": ";
            for (int j = 0; j < *system_->dimX_; j++) {
                clog << etaXs_[i][j] << " ";
            }
            clog << '\n';
        }
        clog << "tau_:\n";
        for (size_t i = 0; i < etaXs_.size(); i++) {
            clog << "abstraction " << i << ": " << *tau_[i] << '\n';
        }
        printVec(Xs_, "X");
        printVec(Os_, "O");
        cout << "U:\n";
        U_->printInfo(1);
        printVec(Gs_, "G");

        checkMakeDir("plotting");
        Xs_[0]->writeToFile("plotting/X.bdd");
        Os_[*system_->numAbs_-1]->writeToFile("plotting/O.bdd");
        Gs_[*system_->numAbs_-1]->writeToFile("plotting/G.bdd");
        checkMakeDir("O");
        saveVec(Os_, "O/O");
        checkMakeDir("G");
        saveVec(Gs_, "G/G");
    }

    /*! Saves a snapshot of a controller C and its domain Z into the sequence of final controllers and controller domains.
        \param[in] ab	0-index of the abstraction which the controller and controller domain that should be saved belong to.
    */
    void saveCZ(int ab) {
        if (Zs_[ab]->symbolicSet_ == ddmgr_->bddZero()) {
            return;
        }

        SymbolicSet* C = new SymbolicSet(*Cs_[ab]);
        C->symbolicSet_ = Cs_[ab]->symbolicSet_;
        finalCs_.push_back(C);
        SymbolicSet* Z = new SymbolicSet(*Xs_[ab]);
        Z->symbolicSet_ = Zs_[ab]->symbolicSet_;
        finalZs_.push_back(Z);
        finalAbs_.push_back(ab);
    }

    /*! Reads transition relation BDDs from file and saves them into Ts_. */
    void loadTs() {
        for (int i = 0; i < *system_->numAbs_; i++) {
            string Str = "T/T";
            Str += std::to_string(i+1);
            Str += ".bdd";
            char Char[20];
            size_t Length = Str.copy(Char, Str.length() + 1);
            Char[Length] = '\0';
            SymbolicSet T(*ddmgr_, Char);
            Ts_[i]->symbolicSet_ = T.symbolicSet_;

            Ts_[i]->printInfo(1);
        }
        clog << "Ts_ read from file.\n";
    }
};
}

#endif /* ADAPTABSREACH_HH_ */
