/*! \file Adaptive.hh
    Contains the Adaptive class.
*/

#ifndef ADAPTIVE_HH_
#define ADAPTIVE_HH_

#include <cstdio>
#include <vector>

#include "SymbolicModelGrowthBound.hh"
#include "FixedPoint.hh"
#include "Helper.hh"
#include "System.hh"


using std::clog;
using std::freopen;
using std::string;
using std::vector;
using namespace helper;

namespace scots {

/*! \class Adaptive
    \brief A class that does abstraction-based synthesis with adaptive time- and state space-gridding.
*/
class Adaptive {
public:

    Cudd* ddmgr_; /*!< A single manager object common to all BDDs used in the program. */
    System* system_; /*!< Contains abstraction parameters. */
    int readXX_; /*!< Whether XXs_ is computed or read from file. */
    int readAbs_; /*!< Whether Abs_ is computed or read from file. */

    int alignment_; /*!< Affects how XX is defined, determined by etaRatio_. */
    vector<double*> etaXs_; /*!< *system_->numAbs_ x *system_->dimX_ matrix of state space grid spacings. */
    vector<double*> tau_; /*!< *system_->numAbs_ x 1 matrix of time steps. */

    vector<SymbolicSet*> Xs_; /*!< The *system_->numAbs_ "pre" state space abstractions, coarsest (0) to finest. */
    vector<SymbolicSet*> Os_; /*!< Instance of *Xs_[i] containing unsafe (obstacle) states. */
    vector<SymbolicSet*> Zs_; /*!< Instance of *Xs_[i] containing winning states. */
    vector<SymbolicSet*> X2s_; /*!< The *system_->numAbs_ "post" state space abstractions, coarsest (0) to finest. */
    vector<SymbolicSet*> XXs_; /*!< The *system_->numAbs_ - 1 mappings between consecutive state space abstractions for which membership implies that the finer cell is a subset of the coarser cell. */
    SymbolicSet* U_; /*!< The single input space abstraction. */
    vector<SymbolicSet*> Cs_; /*!< Controller \subseteq *Xs_[i] x *U_. */

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

    /*!	Constructor for an Adaptive object.
     *  \param[in]  logFile     Filename of program log.
     */
    Adaptive(char* logFile) {
        freopen(logFile, "w", stderr);
        clog << logFile << '\n';
    }

    /*! Destructor for an Adaptive object. */
    ~Adaptive() {
        deleteVecArray(etaXs_);
        deleteVec(tau_);
        deleteVec(Xs_);
        deleteVec(Os_);
        deleteVec(Zs_);
        deleteVec(X2s_);
        deleteVec(XXs_);
        delete U_;
        deleteVec(Cs_);
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
        fclose(stderr);
        delete ddmgr_;        
    }

    /*! Initializes data members.
     *  \param[in]  system      Contains abstraction parameters.
     *  \param[in]	readXX		Whether initializeXXs should be done by construction (0) or reading from files (1).
     *  \param[in]	readAbs		Whether initializeAbs should be done by construction (0) or reading from files (1).
     *  \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
     */
    template<class O_type>
    void initialize(System* system, int readXX, int readAbs, O_type addO) {
        ddmgr_ = new Cudd;
        system_ = system;
        readXX_ = readXX;
        readAbs_ = readAbs;

        initializeAlignment();
        initializeEtaTau();
        initializeSolvers();

        initializeBDDs(addO);
    }

    /*! Initializes BDD-related data members.
     *  \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
     */
    template<class O_type>
    void initializeBDDs(O_type addO) {

        initializeXX2UOZCs(addO);

        TicToc tt;
        tt.tic();
        initializeXXs();
        clog << "------------------------------------------------initializeXXs: ";
        tt.toc();

        initializeNumBDDVars();
        initializePermutes();
        initializeCubes();
        initializeNotVars();

        saveVerify();
    }

    /*! Prints information to the console/log file for the user, and saves some other information. */
    void saveVerify() {
        printEtaX();
        printTau();
        printVec(Xs_, "X");
        printVec(Os_, "O");
        printVec(XXs_, "XX");
        cout << "U:\n";
        U_->printInfo(1);
        checkMakeDir("plotting");
        Xs_[0]->writeToFile("plotting/X.bdd");
        Os_[*system_->numAbs_-1]->writeToFile("plotting/O.bdd");
        checkMakeDir("O");
        saveVec(Os_, "O/O");
    }

//    template<class G_type, class I_type, class O_type, class S_type>
//    void initializeReachAndStay(G_type addG, I_type addI, O_type addO, S_type addS) {
//        if (stage_ != 1) {
//            error("Error: initializeReachAndStay called out of order.\n");
//        }
//        stage_ = 2;

//        TicToc tt;
//        tt.tic();
//        initializeSpec(&Gs_, addG);
//        clog << "Gs_ initialized.\n";
//        initializeSpec(&Is_, addI);
//        clog << "Is_ initialized.\n";
//        initializeSpec(&Os_, addO);
//        clog << "Os_ initialized.\n";
//        initializeSpec(&Ss_, addS);
//        clog << "Ss_ initialized\n";
//        clog << "------------------------------------------------initializeSpec on Gs_, Is_, Os_, Ss_: ";
//        tt.toc();

//        // checking that specification is valid
//        for (int i = 0; i < *system_->numAbs_; i++) {
//            if ((Gs_[i]->symbolicSet_ & Os_[i]->symbolicSet_) != ddmgr_->bddZero()) {
//                error("Error: G and O have nonzero intersection.\n");
//            }
//            if ((Is_[i]->symbolicSet_ & Os_[i]->symbolicSet_) != ddmgr_->bddZero()) {
//                error("Error: I and O have nonzero intersection.\n");
//            }
//            if ((Ss_[i]->symbolicSet_ & Os_[i]->symbolicSet_) != ddmgr_->bddZero()) {
//                error("Error: S and O have nonzero intersection.\n");
//            }
//        }
//        clog << "No obstacle problem with specification.\n";

//        // maximal fixed point starts with whole set
//        for (int i = 0; i < *system_->numAbs_; i++) {
//            Zs_[i]->symbolicSet_ = ddmgr_->bddOne();
//        }
//        clog << "Zs_ modified to BDD-1.\n";

//        for (int i = 0; i < *system_->numAbs_; i++) {
//            SymbolicSet* infZ = new SymbolicSet(*Xs_[i]);
//            infZs_.push_back(infZ);
//        }
//        clog << "infZs_ initialized to empty.\n";
//    }






//    int reachAndStay(int startAbs, int minToGoCoarser, int minToBeValid, int earlyBreak, int verbose = 1) {
//        if (stage_ != 3) {
//            error("Error: reachAndStay called out of order.\n");
//        }
//        TicToc tt;
//        tt.tic();


//    }


//    /*! Maximal fixed-point for an entire abstraction. */
//    void nu(BDD T, int* curAbs) {
//        int iter = 1;
//        while (1) {
//            clog << ".";
//            // get pre of current abtraction's Z disjuncted with projection of converged Z from previous abstraction
//            Cs_[curAbs]->symbolicSet_ = preC(Zs_[curAbs]->symbolicSet_ | infZs_[curAbs]->symbolicSet_, curAbs);
//            // conjunct with safe set (maximal fixed point)
//            Cs_[curAbs]->symbolicSet_ &= Ss_[curAbs]->symbolicSet_;
//            // project onto Xs_[curAbs]
//            BDD Z = Cs_[curAbs]->symbolicSet_.ExistAbstract(*notXvars_[curAbs]);
//        }
//    }



    /*! Calculates the enforceable predecessor of the given set with respect to the transition relation at the specified level of abstraction.
     *  \param[in]  Z           The winning set.
     *  \param[in]  T           The transition relation.
     *  \param[in]  TT          The transition relation with post states existentially abstracted.
     *  \param[in]  curAbs      0-index of the current abstraction.
     *  \return     BDD containing {(x,u)} for which all post states are in Z.
     */
    BDD preC(BDD Z, BDD T, BDD TT, int curAbs) {
        // swap to X2s_[curAbs]
        BDD Z2 = Z.Permute(permutesXtoX2_[curAbs]);
        // posts outside of winning set
        BDD nZ2 = !Z2;
        // {(x,u)} with some posts outside of winning set
        BDD Fbdd = T.AndAbstract(nZ2, *notXUvars_[curAbs]);
        // {(x,u)} with no posts outside of winning set
        BDD nF = !Fbdd;
        // get rid of junk
        BDD preF = TT.AndAbstract(nF, *notXUvars_[curAbs]);

        return preF;
    }

    BDD nPreC(BDD Z, BDD T, BDD TT, int curAbs) {
        // swap to X2s_[curAbs]
        BDD Z2 = Z.Permute(permutesXtoX2_[curAbs]);
        // {(x,u)} with some posts in winning set
        BDD Fbdd = T.AndAbstract(Z2, *notXUvars_[curAbs]);
        // {(x,u)} with no posts in winning set
        BDD nF = !Fbdd & Xs_[curAbs]->symbolicSet_ & U_->symbolicSet_;
//        // get rid of junk
//        BDD preF = TT.AndAbstract(nF, *notXUvars_[curAbs]);
        // project to X
        BDD preZ = nF.ExistAbstract(*notXvars_[curAbs]);
        // second complement
        BDD nPreZ = !preZ & Xs_[curAbs]->symbolicSet_;

        SymbolicSet Z2SS(*X2s_[curAbs]);
        Z2SS.addGridPoints();
        Z2SS.symbolicSet_ &= Z2;
        cout << "Z2:\n";
        Z2SS.printInfo(1);

        SymbolicSet FbddSS(*Cs_[curAbs]);
        FbddSS.symbolicSet_ = Fbdd;
        cout << "Fbdd:\n";
        FbddSS.printInfo(1);

        SymbolicSet nFSS(*Cs_[curAbs]);
        nFSS.symbolicSet_ = nF;
        cout << "nF:\n";
        nFSS.printInfo(1);

        SymbolicSet preZSS(*Xs_[curAbs]);
        preZSS.symbolicSet_ = preZ;
        cout << "preZ:\n";
        preZSS.printInfo(1);

        SymbolicSet nPreZSS(*Xs_[curAbs]);
        nPreZSS.symbolicSet_ = nPreZ;
        cout << "nPreZ:\n";
        nPreZSS.printInfo(1);

        return nPreZ;
    }

    /*! Initializes the BDDs useful for existential abstraction. Must be called only after initializing Xs, U, and X2s. */
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

//    /*!	Inner-approximates a set of states in a coarser abstraction with a finer abstraction.
//        Optimal only if each element of etaRatio_ is 2 (or a greater power of; not 1!).
//        \param[in]      Zc      Winning states in the coarser abstraction.
//        \param[in,out]  Zf      Winning states in the finer abstraction.
//        \param[in]      c       0-index of the coarser abstraction.
//    */
//    void innerFinerAll2(SymbolicSet* Zc, SymbolicSet* Zf, int c) {

//        int dimX = *system_->dimX_;
//        double* etaX = etaXs_[c+1];

//        auto f = [Zc, dimX, etaX](double* x)->bool {
//            std::vector<double> exPlus(x, x + dimX);
//            std::vector<double> exMinus(x, x + dimX);
//            for (int i = 0; i < dimX; i++) {
//                exPlus[i] += etaX[i];
//                exMinus[i] -= etaX[i];

////                clog << exPlus[i] << ' ' << exMinus[i] << '\n';
//            }
////            if (Zc->isElement(exPlus) && Zc->isElement(exMinus)) {
////                clog << "x: " << x[0] << ' ' << x[1] << '\n';
////            }
//            return (Zc->isElement(exPlus) && Zc->isElement(exMinus));

//        };
//        int iter = Zf->addByFunction(f);

//        clog << "iterations: " << iter << '\n';

//    }

    /*!	Inner-approximates a set of states in a finer abstraction with a coarser abstraction.
        This version takes advantage of the fact that cells in the two abstractions are aligned,
        which occurs if all elements of etaRatio_ are a power of 3.
        \param[in,out]  Zc      Winning states in the coarser abstraction.
        \param[in]      Zf      Winning states in the finer abstraction.
        \param[in]      c       0-index of the coarser abstraction.
        \return         1 if Zc grows; 0 otherwise
    */
    int innerCoarserAligned(SymbolicSet* Zc, SymbolicSet* Zf, int c) {
        BDD nQ = !((!(XXs_[c]->symbolicSet_)) | Zf->symbolicSet_);
        BDD cand = (!(nQ.ExistAbstract(*notXvars_[c]))) & Xs_[c]->symbolicSet_;
        if (cand <= Zc->symbolicSet_) {
            return 0;
        }
        else {
            Zc->symbolicSet_ = cand;
            return 1;
        }
    }
//    int innerCoarserAligned(SymbolicSet* Zc, SymbolicSet* Zf, int c) {

//        /* Use numFiner to check if the right number of finer cells are in Zf for a particular corresponding cell in Zc. */
//        int numFiner = 1;
//        for (int i = 0; i < *system_->dimX_; i++) {
//            numFiner *= system_->etaRatio_[i];
//        }

//        SymbolicSet Qcf(*XXs_[c]);
//        Qcf.symbolicSet_ = XXs_[c]->symbolicSet_ & Zf->symbolicSet_;

//        SymbolicSet Qc(*Zc);
//        Qc.symbolicSet_ = Qcf.symbolicSet_.ExistAbstract(*notXvars_[c]); // & S1
//        Qc.symbolicSet_ &= !(Zc->symbolicSet_); /* don't check states that are already in Zc */

//        int* QcMintermWhole;
//        SymbolicSet Ccf(Qcf);
//        BDD result = ddmgr_->bddZero();
//        int* QcMinterm = new int[Qc.nvars_];

//        for (Qc.begin(); !Qc.done(); Qc.next()) { // iterate over all coarse cells with any corresponding finer cells in Zf
//            QcMintermWhole = (int*) Qc.currentMinterm();
//            std::copy(QcMintermWhole + Qc.idBddVars_[0], QcMintermWhole + Qc.idBddVars_[0] + Qc.nvars_, QcMinterm);

//            BDD coarseCell = Qc.mintermToBDD(QcMinterm) & Xs_[c]->symbolicSet_; // a particular coarse cell

//            Ccf.symbolicSet_ = Qcf.symbolicSet_ & coarseCell; // corresponding finer cells to the coarse cell

//            if ((Ccf.symbolicSet_.CountMinterm(Ccf.nvars_)) == (numFiner)) { // if there's a full set of finer cells
//                result |= coarseCell;
//            }
//        }

//        delete[] QcMinterm;

//        if (result == ddmgr_->bddZero()) {
//            return 0;
//        }

//        Zc->symbolicSet_ |= result;
//        return 1;
//    }

    /*! Debugging function. Tests the functions for projecting a set from one state space abtraction to another. */
    template<class G_type>
    void testProjections(G_type addG, int which) {
        int c = 0;
        SymbolicSet Zc(*Xs_[c]);
        SymbolicSet Zf(*Xs_[c+1]);

        if (which == 1) {
            cout << "Testing innerCoarserAligned.\n";
            addG(&Zf);
            TicToc tt;
            tt.tic();

            double p[2] = {2, 2};
            Zc.addPoint(p);

//            innerCoarserRupak(&Zc, &Zf, c);
            innerCoarserAligned(&Zc, &Zf, c);

            tt.toc();
        }
        else {
            cout << "Testing innerFinerAligned.\n";

            addG(&Zc);
            TicToc tt;
            tt.tic();

            innerFinerAligned(&Zc, &Zf, c);
            tt.toc();
        }


        cout << "Zf:\n";
        Zf.printInfo(1);
        cout << "Zc:\n";
        Zc.printInfo(1);

        checkMakeDir("test");
        Zc.writeToFile("test/Zc.bdd");
        Zf.writeToFile("test/Zf.bdd");
    }

    /*! Initializes SymbolicModelGrowthBound objects for each abstraction as well as Ts_, TTs_ for use in the fixed points.
        \param[in]	sysNext		Function pointer to equation that evolves system state.
        \param[in]	radNext		Function pointer to equation that computes growth bound.
    */
    template<class sys_type, class rad_type, class X_type, class U_type>
    void computeAbstractions(sys_type sysNext, rad_type radNext, X_type x, U_type u) {

        for (int i = 0; i < *system_->numAbs_; i++) {
            SymbolicSet* T = new SymbolicSet(*Cs_[i], *X2s_[i]);
            Ts_.push_back(T);
        }
        clog << "Ts_ initialized with empty domain.\n";

        TicToc tt;
        tt.tic();
        if (readAbs_ == 0) {
            for (int i = 0; i < *system_->numAbs_; i++) {
                SymbolicModelGrowthBound<X_type, U_type> Ab(Xs_[i], U_, X2s_[i]);
                Ab.computeTransitionRelation(sysNext, radNext, *solvers_[i]);

                Ts_[i]->symbolicSet_ = Ab.transitionRelation_;
            }
            clog << "Ts_ computed by sampling system behavior.\n";
        }
        else {
            loadTs();
        }

        clog << "------------------------------------------------computeAbstractions: ";
        tt.toc();

        if (readAbs_ == 0) {
            checkMakeDir("T");
            saveVec(Ts_, "T/T");
            clog << "Wrote Ts_ to file.\n";
        }

        for (int i = 0; i < *system_->numAbs_; i++) {
            BDD* TT = new BDD;
            *TT = Ts_[i]->symbolicSet_.ExistAbstract(*cubesX2_[i]);
            TTs_.push_back(TT);
        }
        clog << "TTs_ initialized by ExistAbstracting from Ts_.\n";
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

    /*! Removes a series of sets of states from the pre-states of the series of transition relations.
     *  \param[in] vec      Pointer to vector of SymbolicSets whose BDDs are removed from Ts_ and TTs_.
     */
    template<class vec_type>
    void removeFromTs(vec_type* vec) {
        for (int i = 0; i < *system_->numAbs_; i++) {
            Ts_[i]->symbolicSet_ &= !(vec[0][i]->symbolicSet_);
            *TTs_[i] &= !(vec[0][i]->symbolicSet_);
        }
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

    /*! Initializes the vectors of SymbolicSets Xs_, Os_, X2s_, Us_, Zs_, validZs_, Cs_, validCs_.
     *  \param[in]	addO	Function pointer specifying the points that should be added to the obstacle set.
     */
    template<class O_type>
    void initializeXX2UOZCs(O_type addO) {
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
    }

    /*! Initializes alignment_ by examining etaRatio_. */
    void initializeAlignment() {
        int onlyIntegers = 1;
        int inRange = 1;
        for (int i = 0; i < *system_->dimX_; i++) {
            double etaRatio = system_->etaRatio_[i];
            int etaRatioInt = (int) etaRatio;
            if (!(etaRatio == etaRatioInt)) {
                onlyIntegers = 0;
            }
            if (etaRatio <= 0) {
                inRange = 0;
            }
        }
        if (!onlyIntegers || !inRange) {
            error("Error: unsupported etaRatio, which must be only positive integers.\n");
        }
        alignment_ = 1;
        clog << "Initialized alignment.\n";
    }

    /*! Initializes the abstractions' state space grid parameters and time sampling parameters. */
    void initializeEtaTau() {
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

    /*! Initializes SymbolicSets containing the mappings between consecutive state space abstractions. Reads and writes BDDs from/to file depending on readXXs_. */
    void initializeXXs() {
        for (int i = 0; i < *system_->numAbs_ - 1; i++) {
            SymbolicSet* XX = new SymbolicSet(*Xs_[i], *Xs_[i+1]);

            if (readXX_ == 0) {
                mapAbstractions(Xs_[i], Xs_[i+1], XX);
            }
            else {
                string Str = "XX/XX";
                Str += std::to_string(i+1);
                Str += ".bdd";
                char Char[20];
                size_t Length = Str.copy(Char, Str.length() + 1);
                Char[Length] = '\0';
                SymbolicSet XXSet(*ddmgr_, Char);
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
        for (int i = 0; i < *system_->numAbs_; i++) {
            numBDDVars_ += Xs_[i]->nvars_;
            numBDDVars_ += X2s_[i]->nvars_;
        }
        numBDDVars_ += U_->nvars_;

        clog << "Number of BDD variables: " << numBDDVars_ << '\n';
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
//        SymbolicSet U(*U_);
//        U.symbolicSet_ = *cubeU_;
//        clog << "cubeU: " << '\n';
//        U.printInfo(2);

        clog << "Initialized cubes.\n";
    }

    /*! Initializes the arrays of BDD variable IDs that allow a BDD over an X domain to be projected to the identical BDD over the corresponding X2 domain.
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

    /*! Generates a mapping between two consecutive state space abstractions.
     *  \param[in]      Xc          Coarser state space abstraction.
     *  \param[in]      Xf          Finer state space abstraction.
     *  \param[in,out]  XX    Mapping for a finer cell to the coarser cell that is its superset.
     */
    void mapAbstractions(SymbolicSet* Xc, SymbolicSet* Xf, SymbolicSet* XX) {


        int* XfMinterm;
        double* xPoint = new double[Xc->dim_];
        vector<double> XcPoint (Xc->dim_, 0);
        double* XXPoint = new double[XX->dim_];

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
            for (size_t i = 0; i < Xc->dim_; i++) {
                XcPoint[i] = xPoint[i];
            }
            if (!(Xc->isElement(XcPoint))) {
                continue;
            }
            for (size_t i = 0; i < XX->dim_; i++) {
                XXPoint[i] = xPoint[i % Xc->dim_];
            }
            XX->addPoint(XXPoint);
        }
        cout << '\n';

        delete[] xPoint;
        delete[] XXPoint;
    }


    /*! Prints information regarding the abstractions' grid parameters to the log file. */
    void printEtaX() {
        clog << "etaXs_:\n";
        for (size_t i = 0; i < etaXs_.size(); i++) {
            clog << "abstraction " << i << ": ";
            for (int j = 0; j < *system_->dimX_; j++) {
                clog << etaXs_[i][j] << " ";
            }
            clog << '\n';
        }
    }

    /*! Prints information regarding the abstractions' time sampling parameter to the log file. */
    void printTau() {
        clog << "tau_:\n";
        for (size_t i = 0; i < etaXs_.size(); i++) {
            clog << "abstraction " << i << ": " << *tau_[i] << '\n';
        }
    }

    /*! Debugging function. */
    void checkPermutes() {
        for (int i = 0; i < *system_->numAbs_; i++) {
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

        // for *system_->numAbs_ = 2

        SymbolicSet Xs(*Xs_[0], *Xs_[1]);
        SymbolicSet X2s(*X2s_[0], *X2s_[1]);
        SymbolicSet X(Xs, X2s);
        SymbolicSet all(X, *U_);

        all.symbolicSet_ = *cubesX_[0];
        all.printInfo(2);

//        for (int i = 0; i < *system_->numAbs_; i++) {
//            clog << "X " << i << ":\n";
//            BDD thesevars = ddmgr_->bddOne();
//            for (int j = 0; j < *system_->numAbs_; j++) {
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
