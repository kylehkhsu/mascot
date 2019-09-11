#ifndef BLACKBOXREACH_HH_
#define BLACKBOXREACH_HH_

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
    
    enum ReachResult {CONVERGEDVALID, CONVERGEDINVALID, NOTCONVERGED, INITWINNING};
    
    class BlackBoxReach {
    public:
        Cudd* ddmgr_; /*!< A single manager object common to all BDDs used in the program. */
        System* system_; /*!< Contains abstraction parameters. */
//        vector<Environments*> Es_; /*!< Contains a set of parameterized environments. */
        vector<double*> etaXs_; /*!< *system_->numAbs_ x *system_->dimX_ matrix of state space grid spacings. */
        vector<double*> tau_; /*!< *system_->numAbs_ x 1 matrix of time steps. */
        vector<SymbolicSet*> Xs_; /*!< The *system_->numAbs_ "pre" state space abstractions, coarsest (0) to finest. */
        vector<SymbolicSet*> X0s_; /*!< The overapproximation of the set of initial states. Default is the entire state space outside the obstacle. */
        vector<SymbolicSet*> Os_; /*!< Instance of *Xs_[i] containing unsafe (obstacle) states. */
        vector<SymbolicSet*> Es_; /*!< Vector of size *system_->numAbs_ containing the exclusion zone around the obstacles in different layers. */
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
        vector<BDD*> cubesCoarserTs_; /*!< Helper BDD for projection of transitions to coarser layer. */
        vector<SymbolicSet*> Ts_; /*!< Ts_[i] stores the transition relation of abstraction i. */
        vector<SymbolicSet*> TTs_; /*!< TTs_[i] is Ts_[i] with X2s_[i] existentially abtracted. */
        vector<int*> permutesXtoX2_; /*!< To transform a BDD over X variables to an equivalent one over X2 variables. */
        vector<int*> permutesX2toX_; /*!< To transform a BDD over X2 variables to an equivalent one over X variables. */
        vector<int*> permutesCoarser_; /*!< To project a BDD to the next coarser layer. */
        vector<int*> permutesCoarserTs_; /*!< To project a BDD of transition function to the next coarser layer. */
        vector<int*> permutesFiner_; /*!< To project a BDD to the next finer layer. */
        vector<OdeSolver*> solvers_; /*!< ODE solvers (Runge-Katta approximation) for each abstraction time step. */
        OdeSolver* systemSolver_; /*!< ODE solver (Runge-Kutta approximation) used to solve system trajectory (meant to be finer than the abstract solvers). */
        
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
        int verbose_; /*!< 0(default)=do not print intermediate result, 1=print moderate intermediate result on std I/O, 2=print information including the symbolicset details. */
        
        /*!    Constructor for a BlackBoxReach object.
         *  \param[in]  logFile     Filename of program log.
         */
        BlackBoxReach(const char* logFile, int verbose=0) {
            freopen(logFile, "w", stderr);
            clog << logFile << '\n';
            
            abTime_ = 0;
            synTime_ = 0;
            verbose_ = verbose;
        }

        /*! Copy constructor for a BlackBoxReach object. */
        BlackBoxReach(const BlackBoxReach& other) {
            /* same location as pointed by other */
            ddmgr_=other.ddmgr_;
            system_=other.system_;
            etaXs_=other.etaXs_;
            tau_=other.tau_;
            Xs_=other.Xs_;
            X2s_=other.X2s_;
            U_=other.U_;
            numBDDVars_=other.numBDDVars_;
            cubesX_=other.cubesX_;
            cubesX2_=other.cubesX2_;
            cubeU_=other.cubeU_;
            notXUvars_=other.notXUvars_;
            notXvars_=other.notXvars_;
            cubesCoarser_=other.cubesCoarser_;
            cubesCoarserTs_=other.cubesCoarserTs_;
            Ts_=other.Ts_;
            TTs_=other.TTs_;
            computedDs_=other.computedDs_;
            permutesXtoX2_=other.permutesXtoX2_;
            permutesX2toX_=other.permutesX2toX_;
            permutesCoarser_=other.permutesCoarser_;
            permutesCoarserTs_=other.permutesCoarserTs_;
            permutesFiner_=other.permutesFiner_;
            solvers_=other.solvers_;
            systemSolver_=other.systemSolver_;
            finalAbs_=other.finalAbs_;
            uTs_=other.uTs_;
            m_=other.m_;
            p_=other.p_;
            abTime_=other.abTime_;
            synTime_=other.synTime_;
            verbose_=other.verbose_;
            /* fresh copies of elements related to synthesis */
            for (int i = 0; i < (other.X0s_).size(); i++) {
                SymbolicSet* X0 = new SymbolicSet(*other.X0s_[i]);
                X0s_.push_back(X0);
            }
            for (int i = 0; i < (other.Os_).size(); i++) {
                SymbolicSet* O = new SymbolicSet(*other.Os_[i]);
                Os_.push_back(O);
            }
            for (int i = 0; i < (other.Es_).size(); i++) {
                SymbolicSet* E = new SymbolicSet(*other.Es_[i]);
                Es_.push_back(E);
            }
            for (int i = 0; i < (other.Zs_).size(); i++) {
                SymbolicSet* Z = new SymbolicSet(*other.Zs_[i]);
                Zs_.push_back(Z);
            }
            for (int i = 0; i < (other.Cs_).size(); i++) {
                SymbolicSet* C = new SymbolicSet(*other.Cs_[i]);
                Cs_.push_back(C);
            }
            for (int i = 0; i < (other.Gs_).size(); i++) {
                SymbolicSet* G = new SymbolicSet(*other.Gs_[i]);
                Gs_.push_back(G);
            }
            for (int i = 0; i < (other.validZs_).size(); i++) {
                SymbolicSet* validZ = new SymbolicSet(*other.validZs_[i]);
                validZs_.push_back(validZ);
            }
            for (int i = 0; i < (other.validCs_).size(); i++) {
                SymbolicSet* validC = new SymbolicSet(*other.validCs_[i]);
                validCs_.push_back(validC);
            }
            for (int i = 0; i < (other.finalCs_).size(); i++) {
                SymbolicSet* finalC = new SymbolicSet(*other.finalCs_[i]);
                finalCs_.push_back(finalC);
            }
            for (int i = 0; i < (other.finalZs_).size(); i++) {
                SymbolicSet* finalZ = new SymbolicSet(*other.finalZs_[i]);
                finalZs_.push_back(finalZ);
            }
            for (int i = 0; i < (other.Ds_).size(); i++) {
                SymbolicSet* D = new SymbolicSet(*other.Ds_[i]);
                Ds_.push_back(D);
            }
            for (int i = 0; i < (other.innerDs_).size(); i++) {
                SymbolicSet* innerD = new SymbolicSet(*other.innerDs_[i]);
                innerDs_.push_back(innerD);
            }
//            for (int i = 0; i < (other.computedDs_).size(); i++) {
//                SymbolicSet* computedD = new SymbolicSet(*other.computedDs_[i]);
//                computedDs_.push_back(computedD);
//            }
//            for (int i = 0; i < (other.Ts_).size(); i++) {
//                SymbolicSet* T = new SymbolicSet(*other.Ts_[i]);
//                Ts_.push_back(T);
//            }
//            for (int i = 0; i < (other.TTs_).size(); i++) {
//                SymbolicSet* TT = new SymbolicSet(*other.TTs_[i]);
//                TTs_.push_back(TT);
//            }
        }
        /*! Destructor for a BlackBoxReach object. */
        ~BlackBoxReach() {
            clog << "Abstraction construction time: " << abTime_ << " seconds.\n";
            clog << "Controller synthesis time: " << synTime_ << " seconds.\n";
            clog << "Abs + syn time: " << abTime_ + synTime_ << " seconds.\n";
            deleteVecArray(etaXs_);
            deleteVec(tau_);
            deleteVec(Xs_);
            deleteVec(X0s_);
            deleteVec(Os_);
            deleteVec(Es_);
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
            deleteVec(cubesCoarserTs_);
            deleteVec(Ts_);
            deleteVec(TTs_);
            deleteVecArray(permutesXtoX2_);
            deleteVecArray(permutesX2toX_);
            deleteVecArray(permutesCoarser_);
            deleteVecArray(permutesCoarserTs_);
            deleteVecArray(permutesFiner_);
            deleteVec(solvers_);
            delete systemSolver_;
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
        /*! Only synthesis using the existing abstraction.
         *  IMPORTANT: use only after some abstraction has been computed.
         *  \param[in]  p           Parameter controlling amount of finer layer synthesis attempts.
         */
        void plainReach(int p) {
            m_ = p; // max. iterations for consecutive reachability for non-coarsest layers
            clog << "m: " << m_ << '\n';
            
            // synthesize controller
            int ab = 0;
            eagerReachRecurse(ab);
//            ReachResult result = reach(ab);
//            cout << "ReachResult = " << result << ".\n";
            
//            clog << "\nFinal number of controllers: " << finalCs_.size() << '\n';
//            if (verbose_>0) {
//                cout << "\nFinal number of controllers: " << finalCs_.size() << '\n';
//            }
            
//            checkMakeDir("C");
//            saveVec(finalCs_, "C/C");
//            checkMakeDir("Z");
//            saveVec(finalZs_, "Z/Z");
//            checkMakeDir("T");
//            saveVec(Ts_, "T/T");
//            clog << "Wrote Ts_ to file.\n";
            return;
        }
        /*! Refine some particular region around a point of the state space.
         *  \param[in] xarr     Continuous state around which exploration is needed (can be vector or array)
         *  \param[in] r        Radius of exploration
         */
        template<class X_type, class U_type, class sys_type, class rad_type>
        bool exploreAroundPoint(std::vector<double> x, X_type r, U_type u, sys_type sysNext, rad_type radNext) {
            /* sanity check */
            if(!(x.size()==*system_->dimX_)) {
                std::ostringstream os;
                os << "Error: scots::SymbolicSet::isElement(xarr): x must be of size dim_.";
                throw std::invalid_argument(os.str().c_str());
            }
            for (int i=0; i < *system_->dimX_; i++) {
                if (x[i] > system_->ubX_[i] || x[i] < system_->lbX_[i]) {
                    std::ostringstream os;
                    os << "Error: scots::BlackBoxReach::exploreAround(xarr): xarr is outside the state space.";
                    throw std::invalid_argument(os.str().c_str());
                }
            }
//            /* convert xarr to a vector x (reqd by some of the methods) */
//            std::vector<double> x;
//            for (int i=0; i<*system_->dimX_; i++) {
//                x.push_back(xarr[i]);
//            }
            /* find the current exploration level around x */
            int ab;
            for (int i=*system_->numAbs_-1; i >= 0; i--) {
                if (computedDs_[i]->isElement(x)) {
                    ab = i;
                    break;
                }
            }
            if (ab == *system_->numAbs_-1) { /* return false when already explored upto the finest level */
                return false;
            } else { /* otherwise, the exploration happens in the next finer layer */
                ab = ab + 1;
            }
            /* Set up the polytope that represents the area around x */
            int nofHalfSpaces = 2*(*system_->dimX_);
            double* H = new double[nofHalfSpaces*(*system_->dimX_)];
            for (int i=0; i<nofHalfSpaces*(*system_->dimX_); i++) {
                H[i] = 0;
            }
            int j = 0;
            for (int i=0; i<*system_->dimX_; i++) {
                H[j]=-1;
                j += *system_->dimX_;
                H[j]=1;
                j += *system_->dimX_+1;
            }
            double* h = new double[nofHalfSpaces];
            for (int i=0; i<*system_->dimX_; i++) {
                h[2*i] = -(x[i]-r[i]-Xs_[ab]->z_[i]);
                h[2*i+1] = x[i]+r[i]+Xs_[ab]->z_[i];
            }
            /* explore the portion given by polytope {x' | Hx' <= h} */
            Ds_[ab]->clear();
            Ds_[ab]->addPolytope(nofHalfSpaces,H,h,OUTER);
            X_type y;
            computeAbstraction(ab, sysNext, radNext, y, u);
            delete[] H;
            delete[] h;
            return true;
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
            
//            if (!readAbs) {
//                checkMakeDir("T");
//                saveVec(Ts_, "T/T");
//                clog << "Wrote Ts_ to file.\n";
//            }
            
            // synthesize controller
            int ab = 0;
            eagerReachRecurse(ab);
            
//            clog << "\nFinal number of controllers: " << finalCs_.size() << '\n';
//
//            checkMakeDir("C");
//            saveVec(finalCs_, "C/C");
//            checkMakeDir("Z");
//            saveVec(finalZs_, "Z/Z");
//            checkMakeDir("T");
//            saveVec(Ts_, "T/T");
//            clog << "Wrote Ts_ to file.\n";
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
            if (result == INITWINNING) {
                clog << "result: all initial states are winning\n";
                if (verbose_>0) {
                    cout << "result: all initial states are winning.\n";
                }
            } else if (result == CONVERGEDVALID) {
                clog << "result: converged valid\n";
            } else if (result == CONVERGEDINVALID) {
                clog << "result: converged invalid\n";
            } else {
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
            
            if (result == INITWINNING) {
                /* propagate the winning region to all the finer layers (reqd later to check if the initial states are winning in the finest layer) */
                for (int i=ab; i<*system_->numAbs_-1; i++) {
                    finer(Zs_[i],Zs_[i+1],i);
                    Zs_[i+1]->symbolicSet_ |= validZs_[i+1]->symbolicSet_;
                    validZs_[i+1]->symbolicSet_ = Zs_[i+1]->symbolicSet_;
                }
//                /* also check if the winning domain fully covers Es_ in the finest layer */
//                if (!validZs_[*system_->numAbs_-1]->symbolicSet_ & Es_[*system_->numAbs_-1]->symbolicSet_ == ddmgr_->bddZero())
                    return;
//
//                cout << "All initial states are winning but the exclusion region is not winning yet. Continuing...\n";
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
            
            // begin on-the-fly reachability synthesis
            int ab = 0;
            onTheFlyReachRecurse(ab, sysNext, radNext, x, u);
            
//            clog << "\nFinal number of controllers: " << finalCs_.size() << '\n';
//
//            checkMakeDir("C");
//            saveVec(finalCs_, "C/C");
//            checkMakeDir("Z");
//            saveVec(finalZs_, "Z/Z");
//            checkMakeDir("T");
//            saveVec(Ts_, "T/T");
//            checkMakeDir("D");
//            saveVec(computedDs_, "D/D");
//            clog << "Wrote Ts_ to file.\n";
            
            if (verbose_==2) {
                for (int i = 0; i < finalCs_.size(); i++) {
                    std::cout << "Controller " << i+1 << ":"<< '\n';
                    finalCs_[i]->printInfo();
                    std::cout << "Target " << i+1 << '\n';
                    finalZs_[i]->printInfo();
                }
            } else if (verbose_==1) {
                cout << "Number of controllers: " << finalCs_.size() << ".\n";
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
            
            if (verbose_>0) {
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
            if (result == INITWINNING) {
                clog << "result: all initial states are winning\n";
                if (verbose_>0) {
                    cout << "result: all initial states are winning\n";
                }
            }
            if (result == CONVERGEDVALID) {
                clog << "result: converged valid\n";
                if (verbose_>0)
                    cout << "result: converged valid\n";
            }
            else if (result == CONVERGEDINVALID) {
                clog << "result: converged invalid\n";
                if (verbose_>0)
                    cout << "result: converged invalid\n";
            }
            else {
                clog << "result: not converged\n";
                if (verbose_>0)
                    cout << "result: not converged\n";
            }
            
            if (result != CONVERGEDINVALID) {
                saveCZ(ab);
                validZs_[ab]->symbolicSet_ = Zs_[ab]->symbolicSet_;
                validCs_[ab]->symbolicSet_ = Cs_[ab]->symbolicSet_;
                clog << "saved as snapshot, saved to valids\n";
                if (verbose_>0)
                    cout << "saved as snapshot, saved to valids\n";
            }
            else { // result is CONVERGEDINVALID
                Zs_[ab]->symbolicSet_ = validZs_[ab]->symbolicSet_; // reset this layer's progress
                Cs_[ab]->symbolicSet_ = validCs_[ab]->symbolicSet_;
                clog << "reset to valids\n";
                if (verbose_>0)
                    cout << "reset to valids\n";
            }
            
            if (result == INITWINNING) {
                /* propagate the winning domain up to the finest layer (required for the check in the method isInitWinning) */
                for (int i=ab; i<*system_->numAbs_-1; i++) {
                    finer(Zs_[i], Zs_[i+1], i);
                    Zs_[i+1]->symbolicSet_ |= validZs_[i+1]->symbolicSet_;
                    validZs_[i+1]->symbolicSet_ = Zs_[i+1]->symbolicSet_;
                }
                /* also check if the winning domain fully covers Es_ in the finest layer */
//                if (!validZs_[*system_->numAbs_-1]->symbolicSet_ & Es_[*system_->numAbs_-1]->symbolicSet_ == ddmgr_->bddZero())
                    return;
//
//                cout << "All initial states are winning but the exclusion region is not winning yet. Continuing...\n";
                // debug
//                checkMakeDir("X0");
//                saveVec(X0s_, "X0/X0_");
//                checkMakeDir("vZ");
//                saveVec(validZs_, "vZ/Z");
                // debug end
                
            }
            if (result != NOTCONVERGED) { // ab = 0 always converges
                if (ab == *system_->numAbs_ - 1) {
                    return;
                }
                else { // go finer
                    clog << "Going finer\n";
                    if (verbose_>0)
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
                if (verbose_>0)
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
         * NOTE: unlike original MASCOT, we also compute transitions from the obstacle regions. The avoid part of the specificaiton is addressed during the synthesis.
         *  \param[in]  ab      The layer 0-index.
         *  \param[in]  sysNext     System ODE.
         *  \param[in]  radNext     Growth bound ODE.
         *  \param[in]  x           Dummy state point.
         *  \param[in]  u           Dummy input point.
         */
        template<class sys_type, class rad_type, class X_type, class U_type>
        void computeAbstraction(int ab, sys_type sysNext, rad_type radNext, X_type x, U_type u) {
            BDD D = Ds_[ab]->symbolicSet_ & (!computedDs_[ab]->symbolicSet_); // don't repeat sampling
//            D &= !Os_[ab]->symbolicSet_;
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
                
                // Propagate transitions up to the coarsest layer
                SymbolicSet* Tf = new SymbolicSet(*Ts_[ab]);
                Tf->symbolicSet_ = abstraction.transitionRelation_;
                for (int i=ab; i>0; i--) {
                    SymbolicSet* Tc = new SymbolicSet(*Ts_[i-1]);
                    coarserOuterTs(Tc,Tf,i-1);
                    Ts_[i-1]->symbolicSet_ |= Tc->symbolicSet_;
                    TTs_[i-1]->symbolicSet_ = Ts_[i-1]->symbolicSet_.ExistAbstract(*notXUvars_[i-1]);
                    computedDs_[i-1]->symbolicSet_ |= Ts_[i-1]->symbolicSet_.ExistAbstract(*notXvars_[i-1]);
                }
                // Update the corresponding exploration abstraction
                
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
            if (verbose_>0)
                cout << "abstraction: " << ab << '\n';
            while (1) {
                clog << "iteration: " << i << '\n';
                if (verbose_>0)
                    cout << "iteration: " << i << '\n';
                BDD controllablePreZ1 = controllablePre(Zs_[ab]->symbolicSet_ & !(Es_[ab]->symbolicSet_ | Os_[ab]->symbolicSet_), ab);
                BDD controllablePreZ2 = controllablePre(Zs_[ab]->symbolicSet_ & Es_[ab]->symbolicSet_, ab);
                BDD C = (controllablePreZ1 & !Os_[ab]->symbolicSet_) | (controllablePreZ2 & Es_[ab]->symbolicSet_) | Gs_[ab]->symbolicSet_;
                BDD N = C & (!(Cs_[ab]->symbolicSet_.ExistAbstract(*cubeU_)));
                Cs_[ab]->symbolicSet_ |= N;
                Zs_[ab]->symbolicSet_ = C.ExistAbstract(*notXvars_[ab]) | validZs_[ab]->symbolicSet_;
                
                // debug
//                if (verbose_) {
//                    checkMakeDir("InterimDom");
//                    std::string Str="InterimDom/Z_";
//                    Str+=std::to_string(ab);
//                    Str+="_";
//                    Str+=std::to_string(i);
//                    Str+=".bdd";
//                    char Char[20];
//                    size_t Length = Str.copy(Char, Str.length() + 1);
//                    Char[Length] = '\0';
//                    //                vector<SymbolicSet*> Ztemp = Zs_;
//                    //                for (int i=0; i<ab; i++) {
//                    //                    BDD temp = Zs_[i+1]->symbolicSet_;
//                    //                    finer(Zs_[i],Zs_[i+1],i);
//                    //                    Zs_[i+1]->symbolicSet_ = temp & (!(Zs_[i+1]->symbolicSet_));
//                    //                }
//                    BDD temp = Zs_[ab]->symbolicSet_;
//                    if (ab!=0) {
//                        finer(Zs_[ab-1],Zs_[ab],ab-1);
//                        Zs_[ab]->symbolicSet_ = temp & (!(Zs_[ab]->symbolicSet_));
//                    }
//                    Zs_[ab]->writeToFile(Char);
//                    Zs_[ab]->symbolicSet_ = temp;
//                }
                // debug end
                
                
                BDD iw = X0s_[ab]->symbolicSet_ & (!(Zs_[ab]->symbolicSet_));
                if (iw == ddmgr_->bddZero()) {
                    return INITWINNING;
                }
                if (N == ddmgr_->bddZero() && i != 1) {
//                if (N == ddmgr_->bddZero()) { // by kaushik: "i!=1" doesn't make sense as CONVERGEDINVALID will never be true
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
//            // get rid of the obstacles from pre
//            preF &= !Os_[ab]->symbolicSet_;

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
            // debugging
//            checkMakeDir("testing");
//            Zf->writeToFile("testing/Tf0.bdd");
//            Zc->writeToFile("testing/Tc0.bdd");
        }
        
        /*! Calculates the coarser outer projection of Tf and stores it in Tc.
         *  \param[in]  Tc      Contains result of projection (set of transitions).
         *  \param[in]  Tf      Contains set to project (set of transitions).
         *  \param[in]  c       Layer 0-index of Tc.
         */
        void coarserOuterTs(SymbolicSet* Tc, SymbolicSet* Tf, int c) {
            Tc->symbolicSet_ = (Tf->symbolicSet_.ExistAbstract(*cubesCoarserTs_[c])).Permute(permutesCoarserTs_[c]);
        }
        
        /*! Initializes data members.
         *  \param[in]  system      Contains abstraction parameters.
         *  \param[in]  systemTau   The sampling time used for simulating the concrete system trajectory
         *  \param[in]  systemNSubInt   The number of sub intervals in the RK solver for the simulation of the concrete system trajectory
         */
        void initialize(System* system, double systemTau, int systemNSubInt) {
            ddmgr_ = new Cudd;
            system_ = system;
            
            TicToc timer;
            timer.tic();
            
            // This order should be kept.
            initializeEtaTau();
            initializeSolvers(systemTau,systemNSubInt);
            initializeSymbolicSets();
            initializeNumBDDVars();
            initializePermutes();
            initializeCubes();
            initializeNotVars();
            initializeProjs();
            
            abTime_ += timer.toc();
        }
        
        /*! Initializes SymbolicSet data members. */
        void initializeSymbolicSets() {
            for (int i = 0; i < *system_->numAbs_; i++) {
                SymbolicSet* X = new SymbolicSet(*ddmgr_, *system_->dimX_, system_->lbX_, system_->ubX_, etaXs_[i], tau_[i][0]);
                X->addGridPoints();
                Xs_.push_back(X);
            }
            clog << "Xs_ initialized with full domain.\n";
            
            for (int i = 0; i < *system_->numAbs_; i++) {
                SymbolicSet* X0 = new SymbolicSet(*Xs_[i]);
                X0s_.push_back(X0);
            }
            clog << "X0s_ initialized with empty domain.\n";
            
            for (int i = 0; i < *system_->numAbs_; i++) {
                SymbolicSet* X2 = new SymbolicSet(*Xs_[i], 1);
                X2->addGridPoints();
                X2s_.push_back(X2);
            }
            clog << "X2s_ initialized with full domain.\n";
            
            U_ = new SymbolicSet(*ddmgr_, *system_->dimU_, system_->lbU_, system_->ubU_, system_->etaU_, 0);
            U_->addGridPoints();
            clog << "U_ initialized with full domain.\n";
            
            for (int i = 0; i < *system_->numAbs_; i++) {
                SymbolicSet* O = new SymbolicSet(*Xs_[i]);
                Os_.push_back(O);
            }
//            initializeSpec(&Os_, addO, distance);
            clog << "Os_ initialized with empty domain.\n";
            
            for (int i = 0; i < *system_->numAbs_; i++) {
                SymbolicSet* E = new SymbolicSet(*Xs_[i]);
                Es_.push_back(E);
            }
            clog << "Es_ initialized with empty domain.\n";
            
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
            
            for (int i = 0; i < *system_->numAbs_; i++) {
                SymbolicSet* G = new SymbolicSet(*Xs_[i]);
                Gs_.push_back(G);
            }
//            initializeSpec(&Gs_, addG, distance);
            clog << "Gs_ initialized with empty domain.\n";
            
 
            
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
        
        /*! Initializes the goal and obstacle sets accross all the layers.
         *  \param[in]  HO, HG, HI  The vector containing the normal vectors for obstacles, goal and initial sets respectively
         *  \param[in]  ho, hg, hi  Double vectors
         *  \param[in]  distance    The distance metric used to modify the obstacle and the goal (default to 0).
         */
        template<std::size_t s1, std::size_t s2, std::size_t s3, std::size_t s4, std::size_t s5, std::size_t s6>
        bool initializeSpec(std::vector<std::array<double, s1>> HO, std::vector<std::array<double, s2>> ho, std::vector<std::array<double, s3>> HG, std::vector<std::array<double, s4>> hg, std::vector<std::array<double, s5>> HI, std::vector<std::array<double, s6>> hi, double distance) {
            /* initialize specification in the finest layer */
            /* sanity check: # of obstacles */
            if (HO.size()!=ho.size()) {
                std::ostringstream os;
                os << "Error: the number of obstacles cannot be determined.\n";
                throw std::runtime_error(os.str().c_str());
                return false;
            } else {
                for (int i=0; i<HO.size(); i++) {
                    /* sanity check: # of half spaces should match in HO and ho */
                    if (ho[i].size()*(*system_->dimX_) != HO[i].size()) {
                        std::ostringstream os;
                        os << "Error: the polytope for the obstacle " << i << " is ill defined.\n";
                        throw std::runtime_error(os.str().c_str());
                        return false;
                    } else {
                        /* iniitalize the obstacles */
                        size_t p = ho[i].size();
                        Os_[*system_->numAbs_-1]->addPolytope(p,to_native_array(HO[i]),to_native_array(ho[i]),OUTER);
                    }
                }
            }
            //debug:
//            cout << "\tDone with obstacles.\n";
            /* sanity check: # of goals */
            if (HG.size()!=hg.size()) {
                std::ostringstream os;
                os << "Error: the number of goals cannot be determined.\n";
                throw std::runtime_error(os.str().c_str());
                return false;
            } else {
                for (int i=0; i<HG.size(); i++) {
                    /* sanity check: # of half spaces should match in HG and hg */
                    if (hg[i].size()*(*system_->dimX_) != HG[i].size()) {
                        std::ostringstream os;
                        os << "Error: the polytope for the goal " << i << " is ill defined.\n";
                        throw std::runtime_error(os.str().c_str());
                        return false;
                    } else {
                        /* the goals are deflated by "distance" */
                        for (size_t j=0; j<hg[i].size(); j++) {
                            hg[i][j] = hg[i][j] - distance;
                        }
                        size_t p = hg[i].size();
                        /* first clear any previously stored goal set */
//                        Gs_[*system_->numAbs_-1]->clear();
                        Gs_[*system_->numAbs_-1]->addPolytope(p,to_native_array(HG[i]),to_native_array(hg[i]),INNER);
                    }
                }
            }
            //debug:
//            cout << "\tDone with goal.\n";
            /* sanity check: # of initial states */
            if (HI.size()!=hi.size()) {
                std::ostringstream os;
                os << "Error: the number of initial states cannot be determined.\n";
                throw std::runtime_error(os.str().c_str());
                return false;
            } else {
                for (int i=0; i<HI.size(); i++) {
                    /* sanity check: # of half spaces should match in HI and hi */
                    if (hi[i].size()*(*system_->dimX_) != HI[i].size()) {
                        std::ostringstream os;
                        os << "Error: the polytope for the initial state " << i << " is ill defined.\n";
                        throw std::runtime_error(os.str().c_str());
                        return false;
                    } else {
                        /* the initial states are inflated by "distance" */
                        for (size_t j=0; j<hi[i].size(); j++) {
                            hi[i][j] += distance;
                        }
                        size_t p = hi[i].size();
                        /* first clear any previously stored intial states */
//                        X0s_[*system_->numAbs_-1]->clear();
                        X0s_[*system_->numAbs_-1]->addPolytope(p,to_native_array(HI[i]),to_native_array(hi[i]),OUTER);
                    }
                }
            }
            //debug:
//            cout << "\tDone with initial.\n";
            if (X0s_[*system_->numAbs_-1]->symbolicSet_==ddmgr_->bddZero()) {
                /* ignore this environment in this case */
                cout << "Error: the set of initial states cannot be empty.\n";
                return false;
            }
            /* project the specification up to the coarser layers */
            for (int i = *system_->numAbs_-1; i > 0; i--) {
                coarserInner(Gs_[i-1],Gs_[i],i-1);
                coarserOuter(Os_[i-1],Os_[i],i-1);
                coarserOuter(X0s_[i-1],X0s_[i],i-1);
            }
            /* the exclusion regions are obstacles inflated by "distance" in the coarsest layer \ the obstacles in the respective layer (the subtraction will be done later in this method after the projections) */
            for (int i=0; i<HO.size(); i++) {
                for (size_t j=0; j<ho[i].size(); j++) {
                    ho[i][j] += distance;
                }
                size_t p = ho[i].size();
                Es_[0]->addPolytope(p,to_native_array(HO[i]),to_native_array(ho[i]),OUTER);
            }
            Es_[0]->symbolicSet_ &= !Os_[0]->symbolicSet_;
            /* subtract the obstacles of respective layers from the exclusion regions in the coarsest layer */
            for (int i=1; i<*system_->numAbs_; i++) {
                finer(Es_[i-1],Es_[i],i-1);
                Es_[i]->symbolicSet_ &= !Os_[i]->symbolicSet_;
                /* ignore initial states which are blocked by the obstacles or the exclusion region */
                X0s_[i]->symbolicSet_ &= !(Os_[i]->symbolicSet_ | Es_[i]->symbolicSet_);
            }
            //debug:
//            cout << "\tDone with projection.\n";
            /* checking that specification is valid */
            if ((Gs_[*system_->numAbs_-1]->symbolicSet_ & (!(Os_[*system_->numAbs_-1]->symbolicSet_ | Es_[*system_->numAbs_-1]->symbolicSet_))) == ddmgr_->bddZero()) {
                /* ignore this environment in this case */
                cout << "Error: No goal state is outside the obstacles and the exclusion region in the finest layer.\n";
                return false;
            }
            //debug:
//            cout << "\tChecking intersection of goal with obstacle done.\n";
            if (X0s_[*system_->numAbs_-1]->symbolicSet_==ddmgr_->bddZero()) {
                std::ostringstream os;
                /* ignore this environment in this case */
                cout << "Error: no initial state outside the obstacles and the exclusion in the finest layer.\n";
                return false;
            }
            clog << "No obstacle problem with specification.\n";
            verifySave();
            return true;
        }
        /*! Initializes the goal and obstacle sets accross all the layers.
         *  \param[in]  HO, HG, HI  The vector containing the normal vectors for obstacles, goal and initial sets respectively
         *  \param[in]  ho, hg, hi  Double vectors
         */
        template<std::size_t s1, std::size_t s2, std::size_t s3, std::size_t s4, std::size_t s5, std::size_t s6>
        bool initializeSpec(std::vector<std::array<double, s1>> HO, std::vector<std::array<double, s2>> ho, std::vector<std::array<double, s3>> HG, std::vector<std::array<double, s4>> hg, std::vector<std::array<double, s5>> HI, std::vector<std::array<double, s6>> hi) {
            /* when no distance is passed as argument, initialize the distance as 0 and call initializeSpec */
            double distance = 0.0;
            if (!initializeSpec(HO, ho, HG, hg, HI, hi, distance))
                return false;
            else
                return true;
        }
        /*! Computes only the coarsest transitions.
         *  \param[in]  sysNext     System ODE.
         *  \param[in]  radNext     Growth bound ODE.
         *  \param[in]  x           Dummy state point.
         *  \param[in]  u           Dummy input point. */
        template<class sys_type, class rad_type, class X_type, class U_type>
        void initializeAbstraction(sys_type sysNext, rad_type radNext, X_type x, U_type u) {
            Ds_[0]->addGridPoints();
            computeAbstraction(0, sysNext, radNext, x, u);
        }
        
        /*! Computes only the coarsest transitions along with all the auxiliary exploration transitions.
         *  \param[in]  sysNext     System ODE.
         *  \param[in]  radNext     Growth bound ODE.
         *  \param[in]  x           Dummy state point.
         *  \param[in]  u           Dummy input point. */
        template<class sys_type, class rad_type, class X_type, class U_type>
        void initializeAbstractionWithExplore(sys_type sysNext, rad_type radNext, X_type x, U_type u) {
            // start by synthesizing all uTs_
            Ds_[0]->addGridPoints();
            
            TicToc timer;
            timer.tic();
            computeExplorationAbstractions(sysNext, radNext, x, u);
            abTime_ += timer.toc();
            if (verbose_>0)
                cout << "\nComputation time for exploration abstractions: " << timer.toc() <<"\n";
            
            // Ts_[0] is uTs_[0]
            Ts_[0]->symbolicSet_ = uTs_[0]->symbolicSet_;
            TTs_[0]->symbolicSet_ = Ts_[0]->symbolicSet_.ExistAbstract(*notXUvars_[0]);
            computedDs_[0]->addGridPoints();
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
        void initializeSolvers(double systemTau, int systemNSubInt) {
            for (int i = 0; i < *system_->numAbs_; i++) {
                OdeSolver* solver = new OdeSolver(*system_->dimX_, *system_->nSubInt_, *tau_[i]);
                solvers_.push_back(solver);
            }
            systemSolver_ = new OdeSolver(*system_->dimX_, systemNSubInt, systemTau);
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
            
            // Kaushik
            for (int i=1; i<*system_->numAbs_; i++) {
                BDD* varsCoarserT = new BDD[2*(*system_->dimX_ - ones)];
                int ind = 0;
                for (int j = 0; j < *system_->dimX_; j++) {
                    if (system_->etaRatio_[j] == 2) {
                        varsCoarserT[ind] = ddmgr_->bddVar(Xs_[i]->indBddVars_[j][0]);
                        ind++;
                    }
                }
                for (int j = 0; j < *system_->dimX_; j++) {
                    if (system_->etaRatio_[j] == 2) {
                        varsCoarserT[ind] = ddmgr_->bddVar(X2s_[i]->indBddVars_[j][0]);
                        ind++;
                    }
                }
                
                BDD* cubeCoarserT = new BDD;
                *cubeCoarserT = ddmgr_->bddComputeCube(varsCoarserT, NULL, 2*(*system_->dimX_ - ones));
                cubesCoarserTs_.push_back(cubeCoarserT);
                delete[] varsCoarserT;
            }
            for (int i = 1; i < *system_->numAbs_; i++) {
                int* permuteCoarserT = new int[numBDDVars_];
                for (int ind = 0; ind < numBDDVars_; ind++) {
                    permuteCoarserT[ind] = 0;
                }
                
                for (int dim = 0; dim < *system_->dimX_; dim++) {
                    if (system_->etaRatio_[dim] == 2) {
                        for (size_t projInd = 1; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                            permuteCoarserT[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i-1]->indBddVars_[dim][projInd-1];
                            permuteCoarserT[X2s_[i]->indBddVars_[dim][projInd]] = X2s_[i-1]->indBddVars_[dim][projInd-1];
                        }
                    }
                    else if (system_->etaRatio_[dim] == 1){
                        for (size_t projInd = 0; projInd < Xs_[i]->nofBddVars_[dim]; projInd++) {
                            permuteCoarserT[Xs_[i]->indBddVars_[dim][projInd]] = Xs_[i-1]->indBddVars_[dim][projInd];
                            permuteCoarserT[X2s_[i]->indBddVars_[dim][projInd]] = X2s_[i-1]->indBddVars_[dim][projInd];
                        }
                    }
                }
                for (int dim=0; dim < *system_->dimU_; dim++) {
                    for (size_t ind=0; ind<U_->nofBddVars_[dim]; ind++) {
                        permuteCoarserT[U_->indBddVars_[dim][ind]] = U_->indBddVars_[dim][ind];
                    }
                }
                permutesCoarserTs_.push_back(permuteCoarserT);
                
                clog << "permuteCoarserTs " << i << " to " << i-1 << ": ";
                printArray(permuteCoarserT, numBDDVars_);
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
        
        /*! Sets the transition BDDs */
        /* inputes: the vector of the transition relation BDDs and the computedDs_ BDDs */
        void setTransitions(const std::vector<BDD*> T, const std::vector<BDD*> D) {
            if (T.size()<*system_->numAbs_) {
                error("The number of transition BDDs should match the number of abstraction layers.\n");
            } else if (D.size()<*system_->numAbs_) {
                error("The number of computedDs_ BDDs should match the number of abstraction layers.\n");
            } else {
                for (int i=0; i<*system_->numAbs_; i++) {
                    Ts_[i]->setSymbolicSet(*T[i]);
                    TTs_[i]->symbolicSet_ = Ts_[i]->symbolicSet_.ExistAbstract(*notXUvars_[i]);
                    computedDs_[i]->setSymbolicSet(*D[i]);
                }
            }
        }
        
        /*! Clears sets related to the specification and controller synthesis */
        void clear() {
            for (int i = 0; i < X0s_.size(); i++) {
                X0s_[i]->clear();
            }
            for (int i = 0; i < Os_.size(); i++) {
                Os_[i]->clear();
            }
            for (int i = 0; i < Es_.size(); i++) {
                Es_[i]->clear();
            }
            for (int i = 0; i < Zs_.size(); i++) {
                Zs_[i]->clear();
            }
            for (int i = 0; i < Cs_.size(); i++) {
                Cs_[i]->clear();
            }
            for (int i = 0; i < Gs_.size(); i++) {
                Gs_[i]->clear();
            }
            for (int i = 0; i < validZs_.size(); i++) {
                validZs_[i]->clear();
            }
            for (int i = 0; i < validCs_.size(); i++) {
                validCs_[i]->clear();
            }
            finalCs_.clear();
            finalZs_.clear();
            finalAbs_.clear();
            for (int i = 0; i < Ds_.size(); i++) {
                Ds_[i]->clear();
            }
            for (int i = 0; i < innerDs_.size(); i++) {
                innerDs_[i]->clear();
            }
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
            if (verbose_==2) {
                printVec(Xs_, "X");
                printVec(Os_, "O");
                printVec(Es_, "E");
                cout << "U:\n";
                U_->printInfo(1);
                printVec(Gs_, "G");
            }
//            checkMakeDir("plotting");
//            Xs_[0]->writeToFile("plotting/X.bdd");
//            X0s_[*system_->numAbs_-1]->writeToFile("plotting/X0.bdd");
//            Os_[*system_->numAbs_-1]->writeToFile("plotting/O.bdd");
//            Gs_[*system_->numAbs_-1]->writeToFile("plotting/G.bdd");
//            checkMakeDir("O");
//            saveVec(Os_, "O/O");
//            checkMakeDir("G");
//            saveVec(Gs_, "G/G");
        }
        
        /*! Saves a snapshot of a controller C and its domain Z into the sequence of final controllers and controller domains.
         \param[in] ab    0-index of the abstraction which the controller and controller domain that should be saved belong to.
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
        
        /*! Write final synthesis results and plotting data to files */
        void saveFinalResult() {
            checkMakeDir("C");
            saveVec(finalCs_, "C/C");
            checkMakeDir("Z");
            saveVec(finalZs_, "Z/Z");
            
            checkMakeDir("plotting");
            Xs_[0]->writeToFile("plotting/X.bdd");
            X0s_[*system_->numAbs_-1]->writeToFile("plotting/X0.bdd");
            Os_[*system_->numAbs_-1]->writeToFile("plotting/O.bdd");
            Es_[*system_->numAbs_-1]->writeToFile("plotting/E.bdd");
            Gs_[*system_->numAbs_-1]->writeToFile("plotting/G.bdd");
            checkMakeDir("O");
            saveVec(Os_, "O/O");
            checkMakeDir("E");
            saveVec(Es_, "E/E");
            checkMakeDir("G");
            saveVec(Gs_, "G/G");
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
                TTs_[i]->symbolicSet_ = Ts_[i]->symbolicSet_.ExistAbstract(*notXUvars_[i]);
                
                if (verbose_>0)
                    Ts_[i]->printInfo(1);
            }
            clog << "Ts_ read from file.\n";
        }
        /*! Simulate an abstract controlled trajectory (resolve measurement related non-determinism and initial state non-determinism randomly), and simultaneously compute the shortest distance from the safe set boundary.
         input: obstacles is a vector of the two extreme coordinates of the obstacles which are all assumed to be rectangles. Each element of the vector correspond to one obstacle, whose elements are arranged as: {-lb_x1, ub_x1, -lb_x2, ub_x2, ...} where x1, x2, ... are the state variables, and lb, ub represent the lower and upper bound respectively
         input: trajectory is a sequence of points traversed. trajectory[0] = x0 (to be passed) */
        bool simulateAbs(std::vector<std::vector<double>>& trajectory, std::vector<std::vector<double>> obstacles, double& distance) {
            std::vector<double> x; /* current state */
            std::vector<double> u; /* current control input */
            std::vector<double> xu; /* current state-input pair */
            size_t ab, prevAb; /* current and previous abstraction layer */
            /* the ids of the state space variables and input variables in the controller BDD and transition BDD are 0 to dimX_-1 and dimX_ to dimX_ + dimU_ -1, respectively */
            std::vector<size_t> xind;
            std::vector<size_t> xuind;
            for (size_t i=0; i<*system_->dimX_; i++) {
                xind.push_back(i);
                xuind.push_back(i);
            }
            for (size_t i=0; i<*system_->dimU_; i++) {
                xuind.push_back(*system_->dimX_+i);
            }
//            /* Choose one initial state randomly from the set of initial states */
//            X0s_[0]->getRandomGridPoint(&x);
//            trajectory.push_back(x);
            x = trajectory[0]; /* the initial state */
            if (!(X0s_[0]->isElement(x))) {
                cout << "Error: stots::BlackBoxReach::simulateAbs(trajectory, obstacles, distance): the initial state does not match the design value.";
                return false;
            }
            if (!(finalZs_.back()->isElement(x))) {
                cout << "Error: scots::BlackBoxReach::simulateAbs(trajectory, obstacles, distance): the initial state is outside the winning region.";
                return false;
            }
            /* initialize distance between the trajectory and safe set boundary */
            distance = -1;
            int ab_with_min_dist = 0; /* for debug purpose: the layer id which contribute to the minimum distance */
            /* iterate over all controllers */
            for (int i=finalCs_.size()-1; i>=0; i--) {
                /* find the abstraction layer that corresponds to the present controller (assuming all the abstraction layers' eta are different) */
                /* prevAb correspond to the layer used in finalZs_[i-1] */
                bool flag;
                if (i!=0) {
                    for (size_t j=0; j<*system_->numAbs_; j++) {
                        flag = true;
                        for (size_t k=0; k<*system_->dimX_; k++) {
                            if (finalZs_[i-1]->getEta()[k]!=Xs_[j]->getEta()[k]) {
                                flag = false;
                                break;
                            }
                        }
                        if (flag) {
                            prevAb = j;
                            break;
                        }
                    }
                }
                /* ab correspond to the layer used in finalZs_[i] */
                for (size_t j=0; j<*system_->numAbs_; j++) {
                    flag = true;
                    for (size_t k=0; k<*system_->dimX_; k++) {
                        if (finalZs_[i]->getEta()[k]!=Xs_[j]->getEta()[k]) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        ab = j;
                        break;
                    }
                }
                /* for i=0, the goal is the overall goal Gs; otherwise, the goal is to reach the controller domain of i-1 */
                scots::SymbolicSet goal((i==0) ? *Xs_[*system_->numAbs_-1] : *Xs_[prevAb]);
                if (i==0) {
                    goal.setSymbolicSet(Gs_[*system_->numAbs_-1]->symbolicSet_);
                } else {
                    goal.setSymbolicSet(finalZs_[i-1]->symbolicSet_);
                }
                
                while (1) {
                    /* last point in the trajectory is the current state */
                    x = trajectory.back();
                    /* if the current goal is reached, go to the next controller */
                    if (goal.isElement(x))
                        break;
                    /* pick any random valid control input */
                    if (!getRandomMember(finalCs_[i]->setValuedMap(x,xind),u)) {
                        return false;
                    }
                    /* store the indices of state space and control space */
                    xu.clear();
                    for (size_t j=0; j<*system_->dimX_; j++) {
                        xu.push_back(x[j]);
                    }
                    for (size_t j=0; j<*system_->dimU_; j++) {
                        xu.push_back(u[j]);
                    }
                    /* pick a random successor (meaningful when the abstraction is non-determinstic) */
                    std::vector<double> xx;
                    if (!getRandomMember(Ts_[ab]->setValuedMap(xu,xuind),xx)) {
                        return false;
                    }
                    trajectory.push_back(xx);
                    /* update the distacne */
                    /* obstacles are a collection of rectangles */
                    /* compute the box around xx */
                    std::vector<double> lb1, ub1, lb2, ub2;
                    for (size_t j=0; j<*system_->dimX_; j++) {
                        lb1.push_back(xx[j]-(0.5*etaXs_[ab][j]));
                        ub1.push_back(xx[j]+(0.5*etaXs_[ab][j]));
                    }
                    /* iterate over all obstacles (boxes) */
                    for (size_t j=0; j<obstacles.size(); j++) {
                        /* arrange the inputs in suitable form */
                        for (size_t k=0; k<*system_->dimX_; k++) {
                            lb2.push_back(-obstacles[j][2*k]);
                            ub2.push_back(obstacles[j][2*k+1]);
                        }
                        /* compute the distance */
                        if (distance==-1) { /*first time this loop is entered */
                            distance = computeDistance(lb1, ub1, lb2, ub2);
                            ab_with_min_dist = ab;
                        } else {
                            double new_distance = computeDistance(lb1, ub1, lb2, ub2);
                            if (new_distance<distance) {
                                distance = new_distance;
                                ab_with_min_dist = ab;
                            }
                        }
                        lb2.clear();
                        ub2.clear();
                    }
                } /* End of While loop */
            }/* End of For loop over all controllers in finalCs_ */
            cout << "The distance corresponds to layer " << ab_with_min_dist << "\n";
            // debug
//            writeVecToFile(trajectory,"Figures/abs_traj.txt","clean");
            // debug end
            return true;
        }
        /*! Simulate a concrete controlled trajectory, and simultaneously check if it collides with the obstacles. If yes, which point? (module measurement errors.)
         input: obstacles is a vector of the two extreme coordinates of the obstacles which are all assumed to be rectangles. Each element of the vector correspond to one obstacle, whose elements are arranged as: {-lb_x1, ub_x1, -lb_x2, ub_x2, ...} where x1, x2, ... are the state variables, and lb, ub represent the lower and upper bound respectively
         input: trajectory is a sequence of points traversed. trajectory[0] = x0 (to be passed) */
        template<class sys_type, class X_type, class U_type>
        bool simulateSys(sys_type sys_next,
                         X_type xarr, U_type uarr,
                         std::vector<std::vector<double>>& trajectory,
                         std::vector<std::vector<double>> obstacles,
                         std::vector<std::vector<double>> sys_goal,
                         std::vector<double>& unsafeAt) {
            std::vector<double> x; /* current state */
            std::vector<double> u; /* current control input */
            std::vector<double> xu; /* current state-input pair */
            size_t ab, prevAb; /* current and previous abstraction layer */
            /* the ids of the state space variables and input variables in the controller BDD and transition BDD are 0 to dimX_-1 and dimX_ to dimX_ + dimU_ -1, respectively */
            std::vector<size_t> xind;
            std::vector<size_t> xuind;
            for (size_t i=0; i<*system_->dimX_; i++) {
                xind.push_back(i);
                xuind.push_back(i);
            }
            for (size_t i=0; i<*system_->dimU_; i++) {
                xuind.push_back(*system_->dimX_+i);
            }
//            /* Choose one initial state randomly from the set of initial states */
//            X0s_[0]->getRandomGridPoint(&x);
//            trajectory.push_back(x);
            x = trajectory[0]; /* the initial state */
            if (!(X0s_[0]->isElement(x))) {
                // debug
//                if (!(X0s_[0]->isElement(x))) {
//                    cout << "This is a test.";
//                }
                // debug end
                cout << "Error: stots::BlackBoxReach::simulateSys(trajectory, obstacles, distance): the initial state does not match the design value.\n";
                return false;
            }
            if (!(finalZs_.back()->isElement(x))) {
                cout << "Error: scots::BlackBoxReach::simulateSys(trajectory, obstacles, distance): the initial state is outside the winning region.\n";
                return false;
            }
            /* iterate over all controllers */
            for (int i=finalCs_.size()-1; i>=0; i--) {
                /* find the abstraction layer that corresponds to the present controller (assuming all the abstraction layers' eta are different) */
                /* prevAb correspond to the layer used in finalZs_[i-1] */
                bool flag;
                if (i!=0) {
                    for (size_t j=0; j<*system_->numAbs_; j++) {
                        flag = true;
                        for (size_t k=0; k<*system_->dimX_; k++) {
                            if (finalZs_[i-1]->getEta()[k]!=Xs_[j]->getEta()[k]) {
                                flag = false;
                                break;
                            }
                        }
                        if (flag) {
                            prevAb = j;
                            break;
                        }
                    }
                }
                /* ab correspond to the layer used in finalZs_[i] */
                for (size_t j=0; j<*system_->numAbs_; j++) {
                    flag = true;
                    for (size_t k=0; k<*system_->dimX_; k++) {
                        if (finalZs_[i]->getEta()[k]!=Xs_[j]->getEta()[k]) {
                            flag = false;
                            break;
                        }
                    }
                    if (flag) {
                        ab = j;
                        break;
                    }
                }
                

//                if (i!=0) {
//                    scots::SymbolicSet goal(*Xs_[prevAb]);
//                    goal.setSymbolicSet(finalZs_[i-1]->symbolicSet_);
//                }
                /* arrange the inputs in suitable form */
                std::vector<std::vector<double>> lb,ub;
                for (int j=0; j<obstacles.size(); j++) {
                    std::vector<double> l,u;
                    for (size_t k=0; k<*system_->dimX_; k++) {
                        l.push_back(-obstacles[j][2*k]);
                        u.push_back(obstacles[j][2*k+1]);
                    }
                    lb.push_back(l);
                    ub.push_back(u);
                }
                
                
                while (1) {
                    /* last point in the trajectory is the current state */
                    x = trajectory.back();
                    /* if the current goal is reached, go to the next controller */
                    if (i!=0) {
                        scots::SymbolicSet goal(*Xs_[prevAb]);
                        goal.setSymbolicSet(finalZs_[i-1]->symbolicSet_);
                        if (goal.isElement(x))
                            break;
                    } else {
                        bool goal_reached;
                        /* at i=0, check true satisfaction of the goal, i.e. in terms of sys_goal */
                        for (int l=0; l<sys_goal.size(); l++) {
                            goal_reached = true;
                            for (int m=0; m<*system_->dimX_; m++) {
                                if (x[m]< -sys_goal[l][2*m] ||
                                    x[m]> sys_goal[l][2*m+1]) {
                                    goal_reached = false;
                                    break;
                                }
                            }
                            if (goal_reached)
                                break;
                        }
                        if (goal_reached)
                            break;
                    }
                    
                    /* pick any random valid control input */
                    if (!getRandomMember(finalCs_[i]->setValuedMap(x,xind),u)) {
                        /* no input found means the trajectory left the controller domain; this is also treated as unsafe behavior */
                        unsafeAt = x;
                        return false;
                    }
                    /* find the successor state by simulating the system */
                    /* first convert x and u to arrays */
                    for (size_t k=0; k<*system_->dimX_; k++) {
                        xarr[k] = x[k];
                    }
                    for (size_t k=0; k<*system_->dimU_; k++) {
                        uarr[k] = u[k];
                    }
                    sys_next(xarr,uarr,*systemSolver_);
                    /* store xarr as a vector and append to the trajectory */
                    for (size_t k=0; k<*system_->dimX_; k++) {
                        x[k] = xarr[k];
                    }
                    trajectory.push_back(x);
                    //debug
//                    writeVecToFile(trajectory,"Figures/traj.txt","app");
                    /* collision check */
                    bool collision = true;
                    unsafeAt.clear();
                    /* iterate over all obstacles (boxes) */
                    for (size_t j=0; j<obstacles.size(); j++) {
                        for (size_t k=0; k<*system_->dimX_; k++) {
                            if (x[k]<lb[j][k] || x[k]>ub[j][k]) {
                                collision = false;
                                break;
                            }
                        }
                        if (collision) {
                            unsafeAt = x;
                            return false;
                        }
                        collision = true;
                    }
                }
            }
            return true;
        }
        /* check if the initial states are winning or not */
        bool isInitWinning() {
            BDD iw = X0s_[*system_->numAbs_-1]->symbolicSet_ & (!validZs_[*system_->numAbs_-1]->symbolicSet_);
            if (iw==ddmgr_->bddZero()) {
                return true;
            } else {
                return false;
            }
        }
        /* saves a vector of arrays or vectors in a text file in matlab reachable format */
        template<class T>
        bool writeVecToFile(std::vector<T> vec, const string filename, const string mode="clean") {
            if ((mode.compare("clean") != 0) && (mode.compare("app")) != 0) {
                cout << "Error in BlackBoxReach.hh::saveArray: invalid mode. It has to be either clean or app.\n";
                return false;
            }
//            std::string fstr(mode);
            if (mode.compare("clean") == 0)
                remove(filename.c_str());
            std::ofstream fid(filename, std::ofstream::app);
            for (int i=0; i<vec.size(); i++) {
                for (int j=0; j<vec[i].size(); j++) {
                    fid << vec[i][j] << " ";
                }
                fid << "\n";
            }
            fid.close();
            return true;
        }
    private:
        /* get random element from a vector */
        template<class T>
        bool getRandomMember(std::vector<T> vec, T& member) {
            if (vec.size()==0) {
//                cout << "Error: scots::BlackBoxReach::getRandomMember(vec): vec cannot be empty.\n";
                return false;
            }
            member = vec[rand() % vec.size()];
            return true;
        }
        /* length of a T* object */
        template<class T1, std::size_t SIZE>
        inline T1* to_native_array(const std::array<T1,SIZE> a) {
            T1* b = new T1[a.size()];
            for (int i=0; i<a.size(); i++) {
                b[i] = a[i];
            }
            return b;
        }
        /* compute inf-norm distance between two boxes */
        double computeDistance(std::vector<double> lb1, std::vector<double> ub1, std::vector<double> lb2, std::vector<double> ub2) {
            /* sanity check */
            if (!(lb1.size()==ub1.size() && ub1.size()==lb2.size() && lb2.size()==ub2.size()))
//                cout << "The array dimensions do not match.";
                throw "The array dimensions do not match.";
            
            std::vector<double> diff1, diff2;
            /* first assume that the relative position of the boxes are not side by side */
            /* absolute value of lb1-ub2 */
            for (size_t i=0; i<lb1.size(); i++) {
                diff1.push_back(abs(lb1[i]-ub2[i]));
            }
            /* absolute value of lb2-ub1 */
            for (size_t i=0; i<lb2.size(); i++) {
                diff2.push_back(abs(lb2[i]-ub1[i]));
            }
            /* next, adjust if the boxes are side by side */
            for (size_t i=0; i<*system_->dimX_; i++) {
                if ((lb1[i]>=lb2[i] && lb1[i]<=ub2[i]) ||
                    (ub1[i]>=lb2[i] && ub1[i]<=ub2[i]) ||
                    (lb2[i]>=lb1[i] && lb2[i]<=ub1[i]) ||
                    (ub2[i]>=lb1[i] && ub2[i]<=ub1[i])) {
                    diff1[i] = 0;
                    diff2[i] = 0;
                }
            }
            
            /* pointwise minimum of diff1 and diff2 */
            std::vector<double> diff = diff1;
            for (size_t i=0; i<diff1.size(); i++) {
                diff[i] = std::min(diff1[i],diff2[i]);
            }
            /* the inf norm of diff */
            double distance = 0;
            for (size_t i=0; i<diff.size(); i++) {
//                sum += pow(diff[i],2);
                if (distance<diff[i])
                    distance=diff[i];
            }
            return distance;
        }
};
}

#endif /* BLACKBOXREACH_HH_ */
