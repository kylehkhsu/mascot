/*! \file Adaptive.hh
    Contains the Adaptive class.
*/

#ifndef ADAPTIVE_HH_
#define ADAPTIVE_HH_

#include <cstdio>
#include <vector>

#include "SymbolicModelGrowthBound.hh"
#include "FixedPoint.hh"
#include "RungeKutta4.hh"
#include "Profile.hh"
#include "TicToc.hh"
#include "Helper.hh"


using std::clog;
using std::freopen;
using std::string;
using std::vector;
using namespace helper;

namespace scots {

/*! \class Adaptive
    \brief A class that does abstraction-based synthesis with adaptive time- and state space-gridding.
*/
template<class X_type, class U_type>
class Adaptive {
public:

    Cudd ddmgr_; /*!< A single manager object common to all BDDs used in the program. */
    int dimX_; /*!< Dimensionality of state space. */
    double* lbX_; /*!< Lowermost grid point of state space. */
    double* ubX_; /*!< Uppermost grid point of state space. */
    int dimU_; /*!< Dimensionality of input space. */
    double* lbU_; /*!< Lowermost grid point of input space. */
    double* ubU_; /*!< Uppermost grid point of input space. */
    double* etaU_; /*!< Grid spacing of input space abstraction in each dimension. */
    double* etaRatio_; /*!< Ratio between state space grid spacings of consecutive abstractions. */
    double tauRatio_; /*!< Ratio between time steps of consecutive abstractions. */
    int numAbs_; /*!< Number of abstractions of different granularity. */
    int alignment_; /*!< Affects how XX is defined, determined by etaRatio_. */
    int readXX_; /*!< Whether XXs_ is computed or read from file. */
    int readAbs_; /*!< Whether Abs_ is computed or read from file. */

    vector<double*> etaX_; /*!< numAbs_ x dimX_ matrix of state space grid spacings. */
    vector<double*> tau_; /*!< numAbs_ x 1 matrix of time steps. */

    vector<SymbolicSet*> Xs_; /*!< The numAbs_ "pre" state space abstractions, coarsest (0) to finest. */
    vector<SymbolicSet*> Is_; /*!< Instance of *Xs_[i] containing initial states. */
    vector<SymbolicSet*> Os_; /*!< Instance of *Xs_[i] containing unsafe (obstacle) states. */
    vector<SymbolicSet*> Gs_; /*!< Instance of *Xs_[i] containing goal states. */
    vector<SymbolicSet*> Zs_; /*!< Instance of *Xs_[i] containing winning states. */
    vector<SymbolicSet*> validZs_; /*!< Contains winning states that act as savepoints. */
    vector<SymbolicSet*> X2s_; /*!< The numAbs_ "post" state space abstractions, coarsest (0) to finest. */
    vector<SymbolicSet*> XXs_; /*!< The numAbs_ - 1 mappings between consecutive state space abstractions for which membership implies that the finer cell is a subset of the coarser cell. */
    SymbolicSet* U_; /*!< The single input space abstraction. */
    vector<SymbolicSet*> Cs_; /*!< Controller \subseteq *Xs_[i] x *U_. */
    vector<SymbolicSet*> validCs_; /*!< Controllers that act as savepoints. */

    vector<SymbolicSet*> Ss_; /*!< Instance of *Xs_[i] containing possible safe states. */
    vector<SymbolicSet*> infZs_; /*!< Instance of *Xs_[i] containing projection of convergence of previous maximal fixed points. */

    int numBDDVars_; /*!< Total number of BDD variables used by ddmgr_. */

    vector<BDD*> cubesX_; /*!< *cubesX_[i] is a BDD with a single minterm of all 1s over the domain of *Xs_[i]. */
    vector<BDD*> cubesX2_; /*!< *cubesX2_[i] is a BDD with a single minterm of all 1s over the domain of *X2s_[i]. */
    BDD* cubeU_; /*!< A BDD with a single minterm of all 1s over the domain of the input SymbolicSet. */
    vector<BDD*> notXUvars_; /*!< notXUvars_[i] is a BDD over all ddmgr_ variables whose 1-term is DC for *Xs_[i] and *U_ and 1 otherwise. */
    vector<BDD*> notXvars_; /*!< notXvars_[i] is a BDD over all ddmgr_ variables whose 1-term is DC for *Xs_[i] and 1 otherwise. */

    vector<SymbolicSet*> Ts_; /*!< Ts_[i] stores the transition relation of abstraction i. */
    vector<BDD*> TTs_; /*!< TTs_[i] is Ts_[i] with X2s_[i] existentially abtracted. */

    vector<int*> permutesXtoX2_; /*!< To transform a BDD over X variables to an equivalent one over X2 variables. */
    vector<int*> permutesX2toX_; /*!< To transform a BDD over X2 variables to an equivalent one over X variables. */

    vector<OdeSolver*> solvers_; /*!< ODE solvers (Runge-Katta approximation) for each abstraction time step. */
    vector<SymbolicModelGrowthBound<X_type, U_type>*> Abs_; /*!< Abstractions containing the transition relation \subseteq *Xs_[i] x *U_ x *X2s_[i]. */

    vector<SymbolicSet*> finalCs_; /*!< Sequence of controllers that satisfy the specification. */
    vector<SymbolicSet*> finalZs_; /*!< Sequence of domains of finalCs_. */

    int stage_; /*! Helps ensure user calls methods in correct order. */

    /*!	Constructor for an Adaptive object.
        \param[in]	dimX		Dimensionality of the state space.
        \param[in]	lbX			Lowermost grid point of the state space.
        \param[in]	ubX 		Uppermost grid point of the state space.
        \param[in]	etaX 		Coarsest grid spacing of the state space.
        \param[in]	tau 		Coarsest time step.
        \param[in]	dimU		Dimensionality of the input space.
        \param[in]	lbU 		Lowermost grid point of the input space.
        \param[in]	ubU			Uppermost grid point of the input space.
        \param[in] 	etaU		Grid spacing of the input space.
        \param[in]	etaRatio	Ratio between grid spacings of the input space in consecutive abstractions.
        \param[in]	tauRatio	Ratio between time steps of consecutive abstractions.
        \param[in]  nint        Number of sub-intervals in ODE solving per time step.
        \param[in]	numAbs		Number of abstractions.
        \param[in]	readXX		Whether initializeXXs should be done by construction (0) or reading from files (1).
        \param[in]	readAbs		Whether initializeAbs should be done by construction (0) or reading from files (1).
        \param[in]  logFile     Filename of program log.
    */
    Adaptive(int dimX, double* lbX, double* ubX, double* etaX, double tau,
             int dimU, double* lbU, double* ubU, double* etaU,
             double* etaRatio, double tauRatio, int nint,
             int numAbs, int readXX, int readAbs, char* logFile) {
        freopen(logFile, "w", stderr);
        clog << logFile << '\n';

        dimX_ = dimX;
        lbX_ = lbX;
        ubX_ = ubX;
        dimU_ = dimU;
        lbU_ = lbU;
        ubU_ = ubU;
        etaU_ = etaU;
        etaRatio_ = etaRatio;
        tauRatio_ = tauRatio;
        numAbs_ = numAbs;
        readXX_ = readXX;
        readAbs_ = readAbs;
        stage_ = 1;

        initializeAlignment();
        initializeEtaTau(etaX, tau);
        initializeSolvers(nint);
        initializeXX2UZCs();

        TicToc tt;
        tt.tic();
        initializeXXs();
        clog << "------------------------------------------------initializeXXs: ";
        tt.toc();

        initializeNumBDDVars();
        initializePermutes();
        initializeCubes();
        initializeNotVars();

    }

    /*! Destructor for an Adaptive object. */
    ~Adaptive() {
        deleteVecArray(etaX_);
        deleteVec(tau_);
        deleteVec(Xs_);
        deleteVec(Is_);
        deleteVec(Os_);
        deleteVec(Gs_);
        deleteVec(Zs_);
        deleteVec(validZs_);
        deleteVec(X2s_);
        deleteVec(XXs_);
        delete U_;
        deleteVec(Cs_);
        deleteVec(validCs_);
        deleteVec(cubesX_);
        deleteVec(cubesX2_);
        delete cubeU_;
        deleteVec(notXUvars_);
        deleteVec(notXvars_);
        deleteVec(Ts_);
        deleteVec(TTs_);
        deleteVecArray(permutesXtoX2_);
        deleteVecArray(permutesX2toX_);
        deleteVec(solvers_);
        deleteVec(Abs_);
        deleteVec(finalCs_);
        deleteVec(finalZs_);
        deleteVec(Ss_);
        deleteVec(infZs_);
        fclose(stderr);
    }


    /*! Saves and prints to log file some information related to the safety specification. */
    void saveVerifySafe() {
        printEtaX();
        printTau();
        printVec(Xs_, "X");
        printVec(XXs_, "XX");
        printVec(Ss_, "S");
        clog << "U:\n";
        U_->printInfo(1);
        checkMakeDir("plotting");
        Xs_[0]->writeToFile("plotting/X.bdd");
        checkMakeDir("S");
        saveVec(Ss_,"S/S");
    }

    /*! Saves and prints to log file some information related to the reachability/always-eventually specification. */
    void saveVerifyReach() {
        printEtaX();
        printTau();
        printVec(Xs_, "X");
        printVec(XXs_, "XX");
        printVec(Gs_, "G");
        printVec(Os_, "O");
        printVec(Is_, "I");
        cout << "U:\n";
        U_->printInfo(1);
        checkMakeDir("plotting");
        Xs_[0]->writeToFile("plotting/X.bdd");
        Gs_[numAbs_-1]->writeToFile("plotting/G.bdd");
        Os_[numAbs_-1]->writeToFile("plotting/O.bdd");
        Is_[numAbs_-1]->writeToFile("plotting/I.bdd");
        checkMakeDir("G");
        saveVec(Gs_, "G/G");
        checkMakeDir("O");
        saveVec(Os_, "O/O");
    }



    /*! Initializes objects specific to the following specifications: safe.
        \param[in]	addS	Function pointer specifying the points that should be added to the potential safe set.
    */
    template<class S_type>
    void initializeSafe(S_type addS) {
        if (stage_ != 1) {
            error("Error: initializeSafe called out of order.\n");
        }
        stage_ = 2;

        TicToc tt;
        tt.tic();
        initializeSpec(&Ss_, addS);
        clog << "Ss_ initialized\n";
        clog << "------------------------------------------------initializeSpec on Ss_: ";
        tt.toc();

        // empty obstacles
        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* O = new SymbolicSet(*Xs_[i]);
            Os_.push_back(O);
        }

        // maximal fixed point starts with whole set
        for (int i = 0; i < numAbs_; i++) {
            Zs_[i]->symbolicSet_ = ddmgr_.bddOne();
        }
        clog << "Zs_ modified to BDD-1.\n";

        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* infZ = new SymbolicSet(*Xs_[i]);
            infZs_.push_back(infZ);
        }
        clog << "infZs_ initialized to empty.\n";

        saveVerifySafe();
    }

    /*! Initializes objects specific to the following specifications: always-eventually, reach-while-avoid.
        \param[in]	addG	Function pointer specifying the points that should be added to the goal set.
        \param[in]	addI	Function pointer specifying the points that should be added to the initial set.
        \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
    */
    template<class G_type, class I_type, class O_type>
    void initializeReach(G_type addG, I_type addI, O_type addO) {
        if (stage_ != 1) {
            error("Error: initializeReach called out of order.\n");
        }
        stage_ = 2;

        TicToc tt;
        tt.tic();
        initializeSpec(&Gs_, addG);
        clog << "Gs_ initialized.\n";
        initializeSpec(&Is_, addI);
        clog << "Is_ initialized.\n";
        initializeSpec(&Os_, addO);
        clog << "Os_ initialized.\n";
        clog << "------------------------------------------------initializeSpec on Gs_, Is_, Os_: ";
        tt.toc();

        // checking that specification is valid
        for (int i = 0; i < numAbs_; i++) {
            if ((Gs_[i]->symbolicSet_ & Os_[i]->symbolicSet_) != ddmgr_.bddZero()) {
                error("Error: G and O have nonzero intersection.\n");
            }
            if ((Is_[i]->symbolicSet_ & Os_[i]->symbolicSet_) != ddmgr_.bddZero()) {
                error("Error: I and O have nonzero intersection.\n");
            }
        }
        clog << "No obstacle problem with specification.\n";

        // minimal fixed point starts with null
        for (int i = 0; i < numAbs_; i++) {
            Zs_[i]->symbolicSet_ = ddmgr_.bddZero();
        }
        clog << "Zs_ ensured to be empty.\n";

        saveVerifyReach();
    }

    /*!	Writes, should they exist, a sequence of controller and controller domain BDDs to directories 'C' and 'Z' respectively that satisfy the reachability specification.
        \param[in]  startAbs            0-index of the abstraction to start with.
        \param[in]	minToGoCoarser		Minimum number of growing fixed point iterations needed before an attempt to go to a coarser abstraction.
        \param[in]	minToBeValid		Minimum number of growing fixed point iterations needed before a controller is declared valid.
        \param[in]  earlyBreak          If 1, the synthesis ends as soon as I meets the domain of C.
        \param[in]	verbose				If 1, prints additional information during synthesis to the log file.

        \return     1 if controller(s) satisfying specification is/are synthesized; 0 otherwise.
    */
    int reach(int startAbs, int minToGoCoarser, int minToBeValid, int earlyBreak, int verbose = 1) {
        if (stage_ != 3) {
            error("Error: reach called out of order.\n");
        }

        TicToc tt;
        tt.tic();
        int curAbs = startAbs;

        int iter = 1;
        int justCoarsed = 0;
        int iterCurAbs = 1;
        int reached = 0;
        int stop = 0;

        while (1) {
            mu(minToGoCoarser, minToBeValid, earlyBreak, verbose, &curAbs, &iter, &justCoarsed, &iterCurAbs, &reached, &stop);
            if (stop) {
                break;
            }
        }
        if (reached) {
            clog << "Won.\n";

            checkMakeDir("C");
            saveVec(finalCs_, "C/C");
            checkMakeDir("Z");
            saveVec(finalZs_, "Z/Z");
            clog << "----------------------------------------reach: ";
            tt.toc();
            return 1;
        }
        else {
            clog << "Lost.\n";
            clog << "----------------------------------------reach: ";
            tt.toc();
            return 0;
        }
    }

    /*!	Writes, should they exist, a sequence of controller and controller domain BDDs to directories 'C' and 'Z' respectively that satisfy the Buchi box-diamond (aka always-eventually) specification.
     *  Note: of the n resultant controllers in 'C', the n-1th, .., 1st controllers are to reach, and the nth controller is for when the system state is a goal state.
        \param[in]  startAbs            0-index of the abstraction to start with.
        \param[in]	minToGoCoarser		Minimum number of growing fixed point iterations needed before an attempt to go to a coarser abstraction.
        \param[in]	minToBeValid		Minimum number of growing fixed point iterations needed before a controller is declared valid.
        \param[in]  earlyBreak          If 1, reachability ends as soon as I meets the domain of C.
        \param[in]	verbose				If 1, prints additional information during synthesis to the log file.
        \return     1 if controller(s) satisfying specification is/are synthesized; 0 otherwise.
    */
    int alwaysEventually(int startAbs, int minToGoCoarser, int minToBeValid, int verbose = 1) {
        if (stage_ != 3) {
            error("Error: alwaysEventually called out of order.\n");
        }
        int earlyBreak = 0;

        TicToc tt;
        tt.tic();
        int curAbs = startAbs;

        int nuIter = 1;

        while (1) {
            clog << "------------------------------nu iteration: " << nuIter << "------------------------------\n";

            int muIter = 1;
            int justCoarsed = 0;
            int iterCurAbs = 1;
            int reached = 0;
            int muStop = 0;
            while (1) { // eventually
                mu(minToGoCoarser, minToBeValid, earlyBreak, verbose, &curAbs, &muIter, &justCoarsed, &iterCurAbs, &reached, &muStop);
                if (muStop) {
                    break;
                }
            }

            clog << "\nAlwaysEventually: iteration " << nuIter << "\n";

            // always
            // project to the finest abstraction
            for (int i = curAbs; i < numAbs_ - 1; i++) {
                innerFinerAligned(Zs_[i], Zs_[i+1], i);
            }
            Zs_[numAbs_ - 1]->symbolicSet_ &= Xs_[numAbs_ - 1]->symbolicSet_;
            curAbs = numAbs_ - 1;


            BDD preQ = preC(Zs_[curAbs]->symbolicSet_, curAbs); // preC(X2)
            BDD thing = preQ & Gs_[curAbs]->symbolicSet_; // R \cap preC(X2); converges to controller for goal states
            BDD Q = thing.ExistAbstract(*notXvars_[curAbs]); // converges to goal states that are valid for the AE spec (i.e. can continue infinitely)

            if ((Gs_[curAbs]->symbolicSet_ <= Q) && (Q <= Gs_[curAbs]->symbolicSet_)) {
                clog << "Won.\n";

                SymbolicSet* lastZ = finalZs_.back();
                int abs = -1;
                for (int i = 0; i < numAbs_; i++) {
                    int sameEta = 1;
                    for (int j = 0; j < dimX_; j++) {
                        if (etaX_[i][j] != lastZ->eta_[j]) {
                            sameEta = 0;
                        }
                    }
                    int sameTau = tau_[i][0] == lastZ->tau_;

                    if (sameEta && sameTau) {
                        abs = i;
                        break;
                    }
                }

                clog << "Abstraction of E: " << abs << '\n';

                SymbolicSet E(*Xs_[abs]);
                E.symbolicSet_ = lastZ->symbolicSet_ & !(Gs_[abs]->symbolicSet_);
                E.writeToFile("G/E.bdd");

                SymbolicSet* goalC = new SymbolicSet(*Cs_[curAbs]);
                goalC->symbolicSet_ = thing;
                finalCs_.push_back(goalC);

                SymbolicSet* goalZ = new SymbolicSet(*Zs_[curAbs]);
                goalZ->symbolicSet_ = Q;
                finalZs_.push_back(goalZ);

                clog << finalCs_.size()<< '\n';

                checkMakeDir("C");
                saveVec(finalCs_, "C/C");
                checkMakeDir("Z");
                saveVec(finalZs_, "Z/Z");
                clog << "----------------------------------------alwaysEventually: ";
                tt.toc();

                return 1;
            }
            else {
                clog << "Continuing.\n";


                clog << "Updating Gs.\n";
                Gs_[curAbs]->symbolicSet_ = Q;
                for (int i = curAbs; i > 0; i--) {
                    innerCoarserAligned(Gs_[i-1], Gs_[i], i-1);
                }

                clog << "Resetting Zs, Cs.\n";
                for (int i = 0; i < numAbs_; i++) {
                    Zs_[i]->symbolicSet_ = ddmgr_.bddZero();
                    Cs_[i]->symbolicSet_ = ddmgr_.bddZero();
                }

                size_t j = finalCs_.size();
                for (size_t i = 0; i < j; i++) {
                    delete finalCs_.back();
                    finalCs_.pop_back();
                    delete finalZs_.back();
                    finalZs_.pop_back();
                }

                clog << "finalCs_.size(): " << finalCs_.size() << '\n';

            }
            clog << '\n';
            nuIter += 1;
        }
    }

    /*! Writes, should they exist, controllers of specified abstractions that together satisfy a safety specification. */
    void safe() {
        if (stage_ != 3) {
            error("Error: reach called out of order.\n");
        }

        TicToc tt;
        tt.tic();

        for (int curAbs = 0; curAbs < numAbs_; curAbs++) {
            int iter = 1;
            while (1) {
                clog << ".";
                // get pre of current abtraction's Z disjuncted with projection of converged Z from previous abstraction
                Cs_[curAbs]->symbolicSet_ = preC(Zs_[curAbs]->symbolicSet_ | infZs_[curAbs]->symbolicSet_, curAbs);
                // conjunct with safe set (maximal fixed point)
                Cs_[curAbs]->symbolicSet_ &= Ss_[curAbs]->symbolicSet_;
                // project onto Xs_[curAbs]
                BDD Z = Cs_[curAbs]->symbolicSet_.ExistAbstract(*notXvars_[curAbs]);

                if (Z != Zs_[curAbs]->symbolicSet_) { // not converged
                    Zs_[curAbs]->symbolicSet_ = Z; // update and continue
                }
                else { // converged
                    clog << iter << '\n';
                    Zs_[curAbs]->printInfo(1);
                    break;
                }
                iter += 1;
            }
            if (curAbs != numAbs_ - 1) {
                innerFinerAligned(Zs_[curAbs], infZs_[curAbs+1], curAbs); // obtain projection of converged Z onto next, finer abstraction
                Zs_[curAbs+1]->symbolicSet_ &= !(infZs_[curAbs+1]->symbolicSet_); // take away above from the starting Z of the next abstraction
                Ss_[curAbs+1]->symbolicSet_ &= !(infZs_[curAbs+1]->symbolicSet_); // take away same from the S of the next abstraction

                Ts_[curAbs+1]->symbolicSet_ &= !(infZs_[curAbs+1]->symbolicSet_); // don't consider pre states that already have a controller in a coarser abstraction
                *TTs_[curAbs+1] &= !(infZs_[curAbs+1]->symbolicSet_); // same as above
            }
            Xs_[curAbs]->symbolicSet_ = Zs_[curAbs]->symbolicSet_ | infZs_[curAbs]->symbolicSet_; // for verification purposes
        }

        clog << "----------------------------------------safe: ";
        tt.toc();

        checkMakeDir("Z");
        saveVec(Zs_, "Z/Z");
        checkMakeDir("C");
        saveVec(Cs_, "C/C");

        cout << "Domain of all controllers:\n";
        Xs_[numAbs_-1]->printInfo(1);
        Xs_[numAbs_-1]->writeToFile("plotting/D.bdd");
    }

    /*! Calculates the enforceable predecessor of the given set with respect to the transition relation at the specified level of abstraction.
     *  \param[in]  Z           BDD of the winning set.
     *  \param[in]  curAbs      0-index of the current abstraction.
     *  \return     BDD containing {(x,u)} for which all post states are in Z.
     */
    BDD preC(BDD Z, int curAbs) {
        // swap to X2s_[curAbs]
        BDD Z2 = Z.Permute(permutesXtoX2_[curAbs]);
        // posts outside of winning set
        BDD nZ2 = !Z2;
        // {(x,u)} with posts outside of winning set
        BDD Fbdd = Ts_[curAbs]->symbolicSet_.AndAbstract(nZ2, *notXUvars_[curAbs]);
        // {(x,u)} with no posts outside of winning set
        BDD nF = !Fbdd;
        // get rid of junk
        BDD preF = TTs_[curAbs]->AndAbstract(nF, *notXUvars_[curAbs]);
        return preF;
    }

    /*! One iteration in an adaptive minimal fixed point.
        \param[in]		minToGoCoarser		Minimum number of iterations at an abstraction (not coarsest) before attemping to project to the next coarser abstraction.
        \param[in]		minToBeValid		Minimum number of iterations at an abstraction (not finest) after having gone coarser before a controller is declared valid and saved as a backup.
        \param[in]		earlyBreak			Whether the fixed point should end as soon as a strategy exists for the initial state.
        \param[in]		verbose				Whether messages detailing the steps of the algorithm should be printed to the log file.
        \param[in,out]	curAbs				0-index of the abstraction under consideration for the current iteration.
        \param[in,out]	iter				Counter for the total number of iterations of the minimal fixed point.
        \param[in,out]	justCoarsed			Status of the abstraction under consideration.
        \param[in,out]	iterCurAbs			Counter for the consecutive number of iterations for the current abstraction.
        \param[in,out]	reached				Whether the initial state has been declared winning.
        \param[in,out]	stop				Whether the fixed point should end.
    */
    void mu(int minToGoCoarser, int minToBeValid, int earlyBreak, int verbose, int* curAbs, int* iter, int* justCoarsed, int* iterCurAbs, int* reached, int* stop) {
        clog << "current abstraction: " << *curAbs << '\n';
        clog << "mu iteration: " << *iter << '\n';
        clog << "justCoarsed: " << *justCoarsed << '\n';
        clog << "iterCurAbs: " << *iterCurAbs << '\n';
        clog << "reached: " << *reached << '\n';
        clog << "controllers: " << finalCs_.size() << '\n';

        // get pre(Z)
        BDD preF = preC(Zs_[*curAbs]->symbolicSet_, *curAbs);
        // disjunct with goal (minimal fixed point)
        preF |= Gs_[*curAbs]->symbolicSet_;
        // find new {(x,u)}
        BDD N = preF & (!(Cs_[*curAbs]->symbolicSet_.ExistAbstract(*cubeU_)));
        // add new {(x,u)} to C
        Cs_[*curAbs]->symbolicSet_ |= N;
        // project onto Xs_[*curAbs]
        Zs_[*curAbs]->symbolicSet_ = preF.ExistAbstract(*notXvars_[*curAbs]);

        if (((Zs_[*curAbs]->symbolicSet_ & Is_[*curAbs]->symbolicSet_) != ddmgr_.bddZero()) && (*reached == 0)) {
            *reached = 1;
            if (earlyBreak == 1) {
                saveCZ(*curAbs);
                *stop = 1;
                clog << "\nTotal number of controllers: " << finalCs_.size() << '\n';
                return;
            }
        }

        if (*iter != 1) { // first iteration will always fail since we start with an empty Z
            if (N == ddmgr_.bddZero()) { // if iteration failed to yield new (x,u)
                if (*curAbs == numAbs_ - 1) { // if we're in the finest abstraction
                    if (*reached == 1) {
                        saveCZ(*curAbs);
                    }
                    else {
                    }
                    *stop = 1;
                    clog << "\nTotal number of controllers: " << finalCs_.size() << '\n';
                }
                else {
                    if (verbose) {
                        clog << "No new winning states; going finer.\n";
                    }
                    if (*justCoarsed == 1) { // current controller has not been declared valid
                        if (verbose) {
                            clog << "Current controller has not been declared valid.\n";
                            clog << "Removing last elements of finalCs_ and finalZs_.\n";
                            clog << "Resetting Zs_[" << *curAbs << "] and Cs_[" << *curAbs << "] from valids.\n";
                        }

                        // remove the last saved C and F from finalCs_ and Fs_
                        SymbolicSet* C = finalCs_.back();
                        finalCs_.pop_back();
                        delete(C);
                        SymbolicSet* Z = finalZs_.back();
                        finalZs_.pop_back();
                        delete(Z);
                        *justCoarsed = 0;

                        // reset the current winning states and controller
                        Zs_[*curAbs]->symbolicSet_ = validZs_[*curAbs]->symbolicSet_;
                        Cs_[*curAbs]->symbolicSet_ = validCs_[*curAbs]->symbolicSet_;

                    }
                    else {
                        if (verbose) {
                            clog << "Current controller has been declared valid.\n";
                            clog << "Saving Z and C of abstraction " << *curAbs << " into valids, finals.\n";
                            clog << "Finding inner approximation of Zs_[" << *curAbs << "] and copying into validZs_[" << *curAbs+1 << "].\n";
                        }
                        validZs_[*curAbs]->symbolicSet_ = Zs_[*curAbs]->symbolicSet_;
                        validCs_[*curAbs]->symbolicSet_ = Cs_[*curAbs]->symbolicSet_;

                        saveCZ(*curAbs);
                        innerFinerAligned(Zs_[*curAbs], Zs_[*curAbs+1], *curAbs);
                        validZs_[*curAbs+1]->symbolicSet_ = Zs_[*curAbs+1]->symbolicSet_;
                    }
                    *curAbs += 1;
                    *iterCurAbs = 1;
                }
            }
            else { // if there were new (x,u)
                if ((*justCoarsed == 1) && (*iterCurAbs >= minToBeValid)) {
                    if (verbose) {
                        clog << "Current controller now valid.\n";
                    }
                    *justCoarsed = 0;
                }
                if (*curAbs == 0) {
                    if (verbose) {
                        clog << "More new states, already at coarsest gridding, continuing.\n";
                    }
                    *iterCurAbs += 1;
                }
                else {
                    if (*iterCurAbs >= minToGoCoarser) {
                        if (verbose) {
                            clog << "More new states, minToGoCoarser achieved; try going coarser.\n";
                        }
                        int more = innerCoarserAligned(Zs_[*curAbs-1], Zs_[*curAbs], *curAbs-1);

                        if (more == 0) {
                            if (verbose) {
                                clog << "Projecting to coarser gives no more states in coarser; not changing abstraction.\n";
                            }
                            *iterCurAbs += 1;
                        }
                        else {
                            if (verbose) {
                                clog << "Projecting to coarser gives more states in coarser.\n";
                                clog << "Saving Z and C of abstraction " << *curAbs << " into validZs, validCs, finalZs, finalCs.\n";
                                clog << "Saving Z and C of abstraction " << *curAbs-1 << " into validZs_, validCs.\n";
                            }
                            saveCZ(*curAbs);
                            validZs_[*curAbs]->symbolicSet_ = Zs_[*curAbs]->symbolicSet_;
                            validCs_[*curAbs]->symbolicSet_ = Cs_[*curAbs]->symbolicSet_;
                            *justCoarsed = 1;
                            validZs_[*curAbs-1]->symbolicSet_ = Zs_[*curAbs-1]->symbolicSet_;
                            validCs_[*curAbs-1]->symbolicSet_ = Cs_[*curAbs-1]->symbolicSet_;
                            *curAbs -= 1;
                            *iterCurAbs = 1;
                        }
                    }
                    else {
                        *iterCurAbs += 1;
                    }
                }
            }
        }
        else {
            if (verbose) {
                clog << "First iteration, nothing happens.\n";
            }
        }

        clog << "\n";
        *iter += 1;
    }

    /*! Debugging fucntion. Implementation of basic SCOTS safety (single abstraction) in the Adaptive framework. */
    void safeBasicDebug(int verbose = 0) {
        if (numAbs_ != 1) {
            error("Error: For comparison with SCOTS, numAbs needs to be 1.\n");
        }

        int curAbs = 0;

        TicToc tt;
        tt.tic();

        SymbolicSet C(*Xs_[curAbs], *U_);

        int iter = 1;

        while (1) {
            // get pre(Z)
            C.symbolicSet_ = preC(Zs_[curAbs]->symbolicSet_, curAbs);
            // conjunct with safe set (maximal fixed point)
            C.symbolicSet_ &= Ss_[curAbs]->symbolicSet_;
            // project onto Xs_[curAbs]
            BDD Z = C.symbolicSet_.ExistAbstract(*notXvars_[curAbs]);

            clog << "iter: " << iter << '\n';
            iter += 1;

            if (Z != Zs_[curAbs]->symbolicSet_) {
                Zs_[curAbs]->symbolicSet_ = Z;
            }
            else {
                break;
            }
        }

        tt.toc();

        C.printInfo(1);
        checkMakeDir("debug");
        C.writeToFile("debug/C.bdd");

        SymbolicSet X(*Xs_[curAbs]);
        X.symbolicSet_ = C.symbolicSet_.ExistAbstract(U_->getCube());

        clog << "Domain of controller:\n";
        X.printInfo(1);
    }

    /*! Debugging function. Implementation of basic SCOTS reachability (single abstraction) in the Adaptive framework. */
    template<class G_type, class I_type, class O_type>
    void reachBasicDebug(int startAbs, G_type addG, I_type addI, O_type addO, int debug = 0) {

        int curAbs = startAbs;

        SymbolicSet G(*Xs_[curAbs]);
        addG(&G);
        SymbolicSet O(G);
        addO(&O);
        G.symbolicSet_ &= !(O.symbolicSet_);

//        G.printInfo(1);

        SymbolicSet F(*Xs_[curAbs], *U_);
//        F.symbolicSet_ = G.symbolicSet_ & U_->symbolicSet_;
//        F.printInfo(1);
//        U_->printInfo(1);

        SymbolicSet C(*Xs_[curAbs], *U_);

        BDD Z = ddmgr_.bddZero();

        SymbolicSet domain(C);
        domain.addGridPoints();

        int iter = 1;
        while ((C.symbolicSet_ & Is_[curAbs]->symbolicSet_) == ddmgr_.bddZero()) {

            //             project onto Xs_[curAbs]

            Z = F.symbolicSet_.ExistAbstract(*notXvars_[curAbs]);

            if (iter == debug) {
                SymbolicSet Zss(*Xs_[curAbs]);
                Zss.symbolicSet_ = Z;
                clog << "Z = F.symbolicSet_.ExistAbstract(notXvars)\n";
                Zss.printInfo(1);
            }

            // swap to X2s_[curAbs]
            Z = Z.Permute(permutesXtoX2_[curAbs]);

            SymbolicSet X2(*X2s_[curAbs]);
            X2.symbolicSet_ = Z;
            //            X2.printInfo(1);

            BDD nZ = !Z;

            if (iter == debug) {
                X2.symbolicSet_ = nZ;
                clog << "BDD nZ = !Z\n";
                X2.printInfo(1);
            }

            BDD T = Abs_[curAbs]->transitionRelation_;
            BDD Fbdd = T.AndAbstract(nZ, *notXUvars_[curAbs]);


            SymbolicSet Tss(C, X2);
            Tss.symbolicSet_ = T;
            //            Tss.printInfo(1);

            if (iter == debug) {
                SymbolicSet Fss(*Xs_[curAbs], *U_);
                Fss.symbolicSet_ = Fbdd;
                clog << "BDD Fbdd = T.AndAbstract(nZ, notXUvars)\n";
                Fss.printInfo(1);
            }

            BDD nF = !Fbdd;

            if (iter == debug) {
                SymbolicSet nFss(*Xs_[curAbs], *U_);
                nFss.symbolicSet_ = nF;
                clog << "BDD nF = !Fbdd\n";
                nFss.printInfo(1);
            }

            BDD TT = T.ExistAbstract(*cubesX2_[curAbs]);
            BDD preF = TT.AndAbstract(nF, *notXUvars_[curAbs]);

            if (iter == debug) {
                SymbolicSet preFss(*Xs_[curAbs], *U_);
                preFss.symbolicSet_ = preF;
                clog << "BDD preF = TT.AndAbstract(nF, notXUvars)\n";
                preFss.printInfo(1);
            }

            preF |= G.symbolicSet_;

            BDD N = preF & (!(C.symbolicSet_.ExistAbstract(*cubeU_)));

            if (iter == debug) {
                SymbolicSet Nss(*Xs_[curAbs], *U_);
                Nss.symbolicSet_ = N;
                clog << "BDD N = preF & (!(C.symbolicSet_.ExistAbstract(notXvars)))\n";
                Nss.printInfo(1);
            }

            C.symbolicSet_ |= N;

            F.symbolicSet_ = preF;

            iter++;

        }
        clog << "iterations: " << iter << '\n';
        C.printInfo(1);
        checkMakeDir("debug");
        C.writeToFile("debug/C.bdd");

        // normal SCOTS sanity check
        clog << "------------------------------normal SCOTS-----------------------------------\n";
        FixedPoint fp(Abs_[curAbs]);
        SymbolicSet C2(C);
        C2.symbolicSet_ = fp.reach(G.symbolicSet_, Is_[curAbs]->symbolicSet_, 1, 1);
//        C.symbolicSet_ &= domain.symbolicSet_;
        C2.printInfo(1);

        clog << "Cs are equal: " << ((C.symbolicSet_ <= C2.symbolicSet_) && (C2.symbolicSet_ <= C.symbolicSet_)) << "\n";
//        domain.printInfo(2);
    }

    /*! Saves a snapshot of a controller and its domain into the sequence of final controllers and controller domains.
        \param[in] curAbs	0-index of the abstraction which the controller and controller domain that should be saved belong to.
    */
    void saveCZ(int curAbs) {
        SymbolicSet* C = new SymbolicSet(*Cs_[curAbs]);
        C->symbolicSet_ = Cs_[curAbs]->symbolicSet_;
        finalCs_.push_back(C);
        SymbolicSet* Z = new SymbolicSet(*Xs_[curAbs]);
        Z->symbolicSet_ = Zs_[curAbs]->symbolicSet_;
        finalZs_.push_back(Z);
    }

    /*! Initializes the BDDs useful for existential abstraction. Must be called only after initializing Xs, U, and X2s. */
    void initializeNotVars() {
        for (int i = 0; i < numAbs_; i++) {
            BDD* notXUvars = new BDD;
            BDD* notXvars = new BDD;
            *notXUvars = ddmgr_.bddOne();
            for (int j = 0; j < numAbs_; j++) {
                *notXUvars &= *cubesX2_[j];
                if (i != j) {
                    *notXUvars &= *cubesX_[j];
                }
            }
            *notXvars = *notXUvars & *cubeU_;

            notXUvars_.push_back(notXUvars);
            notXvars_.push_back(notXvars);
        }
    }

    /*!	Inner-approximates a set of states in a coarser abstraction with a finer abstraction.
        This version takes advantage of the fact that cells in the two abstractions are aligned,
        which occurs if all elements of etaRatio_ are a power of 3.
        \param[in]      Zc      Winning states in the coarser abstraction.
        \param[in,out]  Zf      Winning states in the finer abstraction.
        \param[in]      c       0-index of the coarser abstraction.
    */
    void innerFinerAligned(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
        BDD Q = XXs_[c]->symbolicSet_ & Zc->symbolicSet_;
        Zf->symbolicSet_ = Q.ExistAbstract(*notXvars_[c+1]) & Xs_[c+1]->symbolicSet_;
    }

    /*!	Inner-approximates a set of states in a coarser abstraction with a finer abstraction.
        Optimal only if each element of etaRatio_ is 2 (or a greater power of; not 1!).
        \param[in]      Zc      Winning states in the coarser abstraction.
        \param[in,out]  Zf      Winning states in the finer abstraction.
        \param[in]      c       0-index of the coarser abstraction.
    */
    void innerFinerAll2(SymbolicSet* Zc, SymbolicSet* Zf, int c) {

        int dimX = dimX_;
        double* etaX = etaX_[c+1];

        auto f = [Zc, dimX, etaX](double* x)->bool {
            std::vector<double> exPlus(x, x + dimX);
            std::vector<double> exMinus(x, x + dimX);
            for (int i = 0; i < dimX; i++) {
                exPlus[i] += etaX[i];
                exMinus[i] -= etaX[i];

//                clog << exPlus[i] << ' ' << exMinus[i] << '\n';
            }
//            if (Zc->isElement(exPlus) && Zc->isElement(exMinus)) {
//                clog << "x: " << x[0] << ' ' << x[1] << '\n';
//            }
            return (Zc->isElement(exPlus) && Zc->isElement(exMinus));

        };
        int iter = Zf->addByFunction(f);

        clog << "iterations: " << iter << '\n';

    }

    /*!	Inner-approximates a set of states in a finer abstraction with a coarser abstraction.
        This version takes advantage of the fact that cells in the two abstractions are aligned,
        which occurs if all elements of etaRatio_ are a power of 3.
        \param[in,out]  Zc      Winning states in the coarser abstraction.
        \param[in]      Zf      Winning states in the finer abstraction.
        \param[in]      c       0-index of the coarser abstraction.
        \return         1 if Zc grows; 0 otherwise
    */
//    int innerCoarserAligned(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
//        BDD nQ = !((!(XXs_[c]->symbolicSet_)) | Zf->symbolicSet_);
//        BDD Zcandidate = (!(nQ.ExistAbstract(*notXvars_[c]))) & Xs_[c]->symbolicSet_;

//        if (Zcandidate <= Zc->symbolicSet_) {
//            return 0;
//        }
//        else {
//            Zc->symbolicSet_ = Zcandidate;
//            return 1;
//        }
//    }
    int innerCoarserAligned(SymbolicSet* Zc, SymbolicSet* Zf, int c) {

        /* Use numFiner to check if the right number of finer cells are in Zf for a particular corresponding cell in Zc. */
        int numFiner = 1;
        for (int i = 0; i < dimX_; i++) {
            numFiner *= etaRatio_[i];
        }

        SymbolicSet Qcf(*XXs_[c]);
        Qcf.symbolicSet_ = XXs_[c]->symbolicSet_ & Zf->symbolicSet_;

        SymbolicSet Qc(*Zc);
        Qc.symbolicSet_ = Qcf.symbolicSet_.ExistAbstract(*notXvars_[c]); // & S1
        Qc.symbolicSet_ &= !(Zc->symbolicSet_); /* don't check states that are already in Zc */

        int* QcMintermWhole;
        SymbolicSet Ccf(Qcf);
        BDD result = ddmgr_.bddZero();
        int* QcMinterm = new int[Qc.nvars_];

        for (Qc.begin(); !Qc.done(); Qc.next()) { // iterate over all coarse cells with any corresponding finer cells in Zf
            QcMintermWhole = (int*) Qc.currentMinterm();
            std::copy(QcMintermWhole + Qc.idBddVars_[0], QcMintermWhole + Qc.idBddVars_[0] + Qc.nvars_, QcMinterm);

            BDD coarseCell = Qc.mintermToBDD(QcMinterm) & Xs_[c]->symbolicSet_; // a particular coarse cell

            Ccf.symbolicSet_ = Qcf.symbolicSet_ & coarseCell; // corresponding finer cells to the coarse cell

            if ((Ccf.symbolicSet_.CountMinterm(Ccf.nvars_)) == (numFiner)) { // if there's a full set of finer cells
                result |= coarseCell;
            }
        }

        delete[] QcMinterm;

        if (result == ddmgr_.bddZero()) {
            return 0;
        }

        Zc->symbolicSet_ |= result;
        return 1;
    }

    /*! Debugging function. Tests the functions for projecting a set from one state space abtraction to another. */
    template<class G_type>
    void testProjections(G_type addG) {
        int c = 0;
        SymbolicSet Zc(*Xs_[c]);
        SymbolicSet Zf(*Xs_[c+1]);

        int which = 1;

        if (which == 1) {
            addG(&Zf);
            TicToc tt;
            tt.tic();

//            innerCoarserRupak(&Zc, &Zf, c);
//            innerCoarserAligned(&Zc, &Zf, c);

            tt.toc();
        }
        else {
            addG(&Zc);
            TicToc tt;
            tt.tic();

            innerFinerAligned(&Zc, &Zf, c);
            tt.toc();
        }


        clog << "Zf:\n";
        Zf.printInfo(1);
        clog << "Zc:\n";
        Zc.printInfo(1);

        checkMakeDir("test");
        Zc.writeToFile("test/Zc.bdd");
        Zf.writeToFile("test/Zf.bdd");
    }

    /*! Initializes SymbolicModelGrowthBound objects for each abstraction as well as Ts_, TTs_ for use in the fixed points.
        \param[in]	sysNext		Function pointer to equation that evolves system state.
        \param[in]	radNext		Function pointer to equation that computes growth bound.
    */
    template<class sys_type, class rad_type>
    void computeAbstractions(sys_type sysNext, rad_type radNext) {
        if (stage_ != 2) {
            error("Error: computeAbstractions called out of order.\n");
        }
        stage_ = 3;

        TicToc tt;
        tt.tic();
        initializeAbs(sysNext, radNext);
        clog << "------------------------------------------------computeAbstractions: ";
        tt.toc();

        initializeTs();
        if (readAbs_ == 0) {
            checkMakeDir("T");
            saveVec(Ts_, "T/T");
        }

        // removing obstacles from transition relation
        for (int i = 0; i < numAbs_; i++) {
            Ts_[i]->printInfo(1);
            Ts_[i]->symbolicSet_ &= !(Os_[i]->symbolicSet_);
            *TTs_[i] &= !(Os_[i]->symbolicSet_);
        }

        clog << "Os_ removed from Ts_, TTs_.\n";
    }

    /*! Initializes a vector of SymbolicSets that is an instance of Xs_ containing points as specified by addSpec.
     *  \param[in] vec      Vector of SymbolicSets to initialize.
     *  \param[in] addSpec  Function pointer specifying the points to add to each element of vec.
     */
    template<class vec_type, class spec_type>
    void initializeSpec(vec_type* vec, spec_type addSpec) {
        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* Xinstance = new SymbolicSet(*Xs_[i]);
            addSpec(Xinstance);
            vec->push_back(Xinstance);
        }
    }

    /*! Initializes the vectors of SymbolicSets Xs_, X2s_, Us_, Zs_, validZs_, Cs_, validCs_.
     */
    void initializeXX2UZCs() {
        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* X = new SymbolicSet(ddmgr_, dimX_, lbX_, ubX_, etaX_[i], tau_[i][0]);
            X->addGridPoints();
            Xs_.push_back(X);
        }
        clog << "Xs_ initialized with full domain.\n";

        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* X2 = new SymbolicSet(*Xs_[i], 1);
            X2->addGridPoints();
            X2s_.push_back(X2);
        }
        clog << "X2s_ initialized with full domain.\n";

        U_ = new SymbolicSet(ddmgr_, dimU_, lbU_, ubU_, etaU_, 0);
        U_->addGridPoints();
        clog << "U_ initialized with full domain.\n";

        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* Z = new SymbolicSet(*Xs_[i]);
            Zs_.push_back(Z);
            SymbolicSet* validZ = new SymbolicSet(*Xs_[i]);
            validZs_.push_back(validZ);
        }
        clog << "Zs_, validZs_ initialized with empty domain.\n";

        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* C = new SymbolicSet(*Xs_[i], *U_);
            Cs_.push_back(C);
            SymbolicSet* validC = new SymbolicSet(*C);
            validCs_.push_back(validC);
        }
        clog << "Cs_, validCs_ initialized with empty domain.\n";
    }

    /*! Initializes Ts_ and TTs_.*/
    void initializeTs() {
        for (int i = 0; i < numAbs_; i++) {
            SymbolicSet* T = new SymbolicSet(*Cs_[i], *X2s_[i]);
            T->symbolicSet_ = Abs_[i]->transitionRelation_;

            BDD* TT = new BDD;
            *TT = T->symbolicSet_.ExistAbstract(*cubesX2_[i]);

            Ts_.push_back(T);
            TTs_.push_back(TT);
        }
    }

    /*! Initializes alignment_ by examining etaRatio_. */
    void initializeAlignment() {
        int isAligned = 1;
        int is2 = 1;
        for (int i = 0; i < dimX_; i++) {
            if (etaRatio_[i] != 2) {
                is2 = 0;
            }
            if ((etaRatio_[i] != 1) && (etaRatio_[i] != 3)) {
                isAligned = 0;
            }
        }
        if (!isAligned && !is2) {
            error("Error: unsupported etaRatio, which must be either a combination of 1s and 3s, or only 2s.\n");
        }
        if (isAligned) {
            alignment_ = 13;
        }
        else if (is2) {
            error("Error: unsupported etaRatio, which must be a combination of 1s and 3s.\n");
            alignment_ = 2;
        }
    }

    /*! Initializes the abstractions' state space grid parameters and time sampling parameters. */
    void initializeEtaTau(double* etaX, double tau) {
        double* etaCur = new double[dimX_];
        for (int i = 0; i < dimX_; i++) {
            etaCur[i] = etaX[i];
        }
        double* tauCur = new double;
        *tauCur = tau;

        for (int i = 0; i < numAbs_; i++) {
            double* etai = new double[dimX_];
            double* taui = new double;
            for (int j = 0; j < dimX_; j++) {
                etai[j] = etaCur[j];
            }
            *taui = *tauCur;
            etaX_.push_back(etai);
            tau_.push_back(taui);

            for (int j = 0; j < dimX_; j++) {
                etaCur[j] /= etaRatio_[j];
            }
            *tauCur /= tauRatio_;
        }
        delete[] etaCur;
        delete tauCur;
    }

    /*! Initializes the Runge-Katta ODE solvers.
     *  \param[in]  nint        Number of desired sub-steps for each time step in the numerical approximation process.
     */
    void initializeSolvers(int nint) {
        for (int i = 0; i < numAbs_; i++) {
            OdeSolver* solver = new OdeSolver(dimX_, nint, *tau_[i]);
            solvers_.push_back(solver);
        }
    }

    /*! Initializes SymbolicSets containing the mappings between consecutive state space abstractions. Reads and writes BDDs from/to file depending on readXXs_. */
    void initializeXXs() {
        for (int i = 0; i < numAbs_ - 1; i++) {
            SymbolicSet* XX = new SymbolicSet(*Xs_[i], *Xs_[i+1]);

            if (readXX_ == 0) {
                mapAbstractions(Xs_[i], Xs_[i+1], XX, i);
            }
            else {
                string Str = "XX/XX";
                Str += std::to_string(i+1);
                Str += ".bdd";
                char Char[20];
                size_t Length = Str.copy(Char, Str.length() + 1);
                Char[Length] = '\0';
                SymbolicSet XXSet(ddmgr_, Char);
                XX->symbolicSet_ = XXSet.symbolicSet_;
            }
            XXs_.push_back(XX);
        }

        if (readXX_ == 0) {
            checkMakeDir("XX");
            saveVec(XXs_, "XX/XX");
        }
    }

    /*! Initializes the number of total distinct BDD variables in use. */
    void initializeNumBDDVars() {
        numBDDVars_ = 0;
        for (int i = 0; i < numAbs_; i++) {
            numBDDVars_ += Xs_[i]->nvars_;
            numBDDVars_ += X2s_[i]->nvars_;
        }
        numBDDVars_ += U_->nvars_;

        clog << "Number of BDD variables: " << numBDDVars_ << '\n';
    }

    /*! Initializes the BDDs that are the precursors to notXvars_ and notXUvars_ */
    void initializeCubes() {
        for (int i = 0; i < numAbs_; i++) {
            BDD* cubeX = new BDD;
            BDD* cubeX2 = new BDD;

            BDD* varsX = new BDD[Xs_[i]->nvars_];
            BDD* varsX2 = new BDD[X2s_[i]->nvars_];

            for (size_t j = 0; j < Xs_[i]->nvars_; j++) {
                varsX[j] = ddmgr_.bddVar(Xs_[i]->idBddVars_[j]);
            }

            for (size_t j = 0; j < X2s_[i]->nvars_; j++) {
                varsX2[j] = ddmgr_.bddVar(X2s_[i]->idBddVars_[j]);
            }

            *cubeX = ddmgr_.bddComputeCube(varsX, NULL, Xs_[i]->nvars_);
            *cubeX2 = ddmgr_.bddComputeCube(varsX2, NULL, X2s_[i]->nvars_);


            cubesX_.push_back(cubeX);
            cubesX2_.push_back(cubeX2);

            delete[] varsX;
            delete[] varsX2;
        }

        cubeU_ = new BDD;
        *cubeU_ = U_->getCube();
//        SymbolicSet U(*U_);
//        U.symbolicSet_ = *cubeU_;
//        clog << "cubeU: " << '\n';
//        U.printInfo(2);
    }

    /*! Initializes the arrays of BDD variable IDs that allow a BDD over an X domain to be projected to the identical BDD over the corresponding X2 domain.
     */
    void initializePermutes() {
        for (int i = 0; i < numAbs_; i++) {
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
    }

    /*! Initializes the transition relations for the abstractions.
     *  \param[in]  sysNext     Specification of system dynamics.
     *  \param[in]  radNext     Specification of state uncertainty growth bound.
     */
    template<class sys_type, class rad_type>
    void initializeAbs(sys_type sysNext, rad_type radNext) {

        for (int i = 0; i < numAbs_; i++) {
            SymbolicModelGrowthBound<X_type, U_type>* Ab = new SymbolicModelGrowthBound<X_type, U_type>(Xs_[i], U_, X2s_[i]);
            if (readAbs_ == 0) {
                Ab->computeTransitionRelation(sysNext, radNext, *solvers_[i]);

            }
            else {
                string Str = "T/T";
                Str += std::to_string(i+1);
                Str += ".bdd";
                char Char[20];
                size_t Length = Str.copy(Char, Str.length() + 1);
                Char[Length] = '\0'; 
                SymbolicSet T(ddmgr_, Char);
                Ab->transitionRelation_ = T.symbolicSet_;
            }

            std::clog << "Number of elements in the transition relation: " << Ab->getSize() << std::endl;

            Abs_.push_back(Ab);
        }
    }

    /*! Generates a mapping between two consecutive state space abstractions.
     *  \param[in]      Xc          Coarser state space abstraction.
     *  \param[in]      Xf          Finer state space abstraction.
     *  \param[in,out]  XX    Mapping for a finer cell to the coarser cell that is its superset.
     *  \param[in]      c           0-index of the coarser abstraction (currently unnecessary/unused).
     */
    void mapAbstractions(SymbolicSet* Xc, SymbolicSet* Xf, SymbolicSet* XX, int c) {
        if (alignment_ == 13) {
            int* XfMinterm;
            double* xPoint = new double[dimX_];
            vector<double> XcPoint (dimX_, 0);
            double* XXPoint = new double[dimX_ + dimX_];

            int totalIter = Xf->symbolicSet_.CountMinterm(Xf->nvars_);
            int iter = 0;
//            int progress = 0;

            cout << "XXs totalIter: " << totalIter << '\n';

            for (Xf->begin(); !Xf->done(); Xf->next()) {
                iter++;
//                cout << iter << '\n';
                if (iter % 50000 == 0) {
                    cout << iter << "\n";
//                    progress++;
                }

                XfMinterm = (int*)Xf->currentMinterm();
                Xf->mintermToElement(XfMinterm, xPoint);
                for (int i = 0; i < dimX_; i++) {
                    XcPoint[i] = xPoint[i];
                }
                if (!(Xc->isElement(XcPoint))) {
                    continue;
                }
                for (int i = 0; i < dimX_+dimX_; i++) {
                    XXPoint[i] = xPoint[i % dimX_];
                }
                XX->addPoint(XXPoint);
            }
            cout << '\n';
            (void)c;

            delete[] xPoint;
            delete[] XXPoint;

        }
//        else if (alignment_ == 2) {
//            int* XfMinterm;
//            double xPoint[dimX_] = {0};
//            double xPointPlus[dimX_] = {0};
//            double xPointMinus[dimX_] = {0};
//            for (Xf->begin(); !Xf->done(); Xf->next()) {
//                XfMinterm = (int*)Xf->currentMinterm();
//                Xf->mintermToElement(XfMinterm, xPoint);
//                for (int i = 0; i < dimX_; i++) {
//                    xPointPlus[i] = xPoint[i] + etaX_[c][i];
//                    xPointMinus[i] = xPoint[i] - etaX_[c][i];
//                }


//            }
//        }
    }

    /*! Throws a logic error along with specified message and closes the log file.
     *  \param[in]  msg     Error message to log to file.
     */
    template<class msg_type>
    void error(msg_type msg) {
        std::ostringstream os;
        os << msg;
        throw std::logic_error(os.str().c_str());
        fclose(stderr);
    }


    /*! Prints information regarding the abstractions' grid parameters to the log file. */
    void printEtaX() {
        clog << "etaX_:\n";
        for (size_t i = 0; i < etaX_.size(); i++) {
            clog << "abstraction " << i << ": ";
            for (int j = 0; j < dimX_; j++) {
                clog << etaX_[i][j] << " ";
            }
            clog << '\n';
        }
    }

    /*! Prints information regarding the abstractions' time sampling parameter to the log file. */
    void printTau() {
        clog << "tau_:\n";
        for (size_t i = 0; i < etaX_.size(); i++) {
            clog << "abstraction " << i << ": " << *tau_[i] << '\n';
        }
    }

    /*! Debugging function. */
    void checkPermutes() {
        for (int i = 0; i < numAbs_; i++) {
            clog << "X:\n";
            Xs_[i]->printInfo(1);
            clog << "X2:\n";
            X2s_[i]->printInfo(1);

            clog << "XtoX2:\n";
            printArray(permutesXtoX2_[i], numBDDVars_);

            clog << "X2toX:\n";
            printArray(permutesX2toX_[i], numBDDVars_);

        }
    }

    /*! Debugging function. */
    void checkCubes() {

        // for numAbs_ = 2

        SymbolicSet Xs(*Xs_[0], *Xs_[1]);
        SymbolicSet X2s(*X2s_[0], *X2s_[1]);
        SymbolicSet X(Xs, X2s);
        SymbolicSet all(X, *U_);

        all.symbolicSet_ = *cubesX_[0];
        all.printInfo(2);

//        for (int i = 0; i < numAbs_; i++) {
//            clog << "X " << i << ":\n";
//            BDD thesevars = ddmgr_.bddOne();
//            for (int j = 0; j < numAbs_; j++) {
//                thesevars &= *cubesX2_[j];
//                if (i != j) {
//                    thesevars &= *cubesX_[j];
//                }
//            }
//            thesevars &= *cubeU_;
//            all.symbolicSet_ = thesevars;
//            all.printInfo(2);
//            clog << "X2 " << i << ":\n";
//            all.symbolicSet_ = *cubesX2_[i];
//            all.printInfo(2);
//        }
    }


    /*! Debugging function. */
    void debug() {
        clog << "\n--------------------debug--------------------\n";
        checkPermutes();
        checkCubes();
    }

};
}

#endif /* ADAPTIVE_HH_ */
