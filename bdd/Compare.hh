/*! \file Compare.hh
    Contains the Compare class.
*/

#ifndef COMPARE_HH_
#define COMPARE_HH_

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

/*! \class Compare
    \brief A class that provides comparison using basic SCOTS.
*/
template<class X_type, class U_type>
class Compare {
public:

    Cudd ddmgr_; /*!< A single manager object common to all BDDs used in the program. */
    int readAb_; /*!< Whether transition relation is computed or read from file. */

    SymbolicSet* X_; /*!< State space. */
    SymbolicSet* G_; /*!< Instance of *X_ containing goal states. */
    SymbolicSet* I_; /*!< Instance of *X_ containing initial states. */
    SymbolicSet* O_; /*!< Instance of *X_ containing unsafe (obstacle) states. */
    SymbolicSet* X2_; /*!< The "post" state space. */
    SymbolicSet* U_; /*!< The input space abstraction. */
    SymbolicSet* C_; /*!< Controller. */
    SymbolicSet* T_; /*!< Transition relation. */

    SymbolicSet* S_; /*!< Instance of *X_ containing possible safe states. */
    SymbolicSet* Z_;

    OdeSolver* solver_; /*!< ODE solver (Runge-Katta approximation) for an abstraction time step. */
    SymbolicModelGrowthBound<X_type, U_type>* Ab_; /*!< Abstraction containing the transition relation. */

    int stage_; /*! Helps ensure user calls methods in correct order. */

    /*!	Constructor for a Compare object.
        \param[in]	dimX		Dimensionality of the state space.
        \param[in]	lbX			Lowermost grid point of the state space.
        \param[in]	ubX 		Uppermost grid point of the state space.
        \param[in]	etaX 		Grid spacing of the state space.
        \param[in]	tau 		Time step.
        \param[in]	dimU		Dimensionality of the input space.
        \param[in]	lbU 		Lowermost grid point of the input space.
        \param[in]	ubU			Uppermost grid point of the input space.
        \param[in] 	etaU		Grid spacing of the input space.
        \param[in]	etaRatio	(Adaptive.hh) Ratio between grid spacings of the input space in consecutive abstractions.
        \param[in]	tauRatio	(Adaptive.hh) Ratio between time steps of consecutive abstractions.
        \param[in]  nint        Number of sub-intervals in ODE solving per time step.
        \param[in]	numAbs		(Adaptive.hh) Number of abstractions.
        \param[in]	readAb		Whether transition relation should constructed (0) or read from file (1).
        \param[in]  logFile     Filename of program log.
    */
    Compare(int dimX, double* lbX, double* ubX, double* etaX, double tau,
             int dimU, double* lbU, double* ubU, double* etaU,
             double* etaRatio, double tauRatio,
             int nint, int numAbs, int readAb, char* logFile) {
        freopen(logFile, "w", stderr);
        clog << logFile << '\n';

        for (int i = 0; i < numAbs - 1; i++) {
            for (int j = 0; j < dimX; j++) {
                etaX[j] = etaX[j] / etaRatio[j];
            }
            tau = tau / tauRatio;
        }

        printEtaX(etaX, dimX);
        printTau(tau);

        readAb_ = readAb;

        X_ = new SymbolicSet(ddmgr_, dimX, lbX, ubX, etaX, tau);
        X_->addGridPoints();
        G_ = new SymbolicSet(*X_);
        I_ = new SymbolicSet(*X_);
        O_ = new SymbolicSet(*X_);
        X2_ = new SymbolicSet(*X_, 1);
        U_ = new SymbolicSet(ddmgr_, dimU, lbU, ubU, etaU, 0);
        U_->addGridPoints();
        C_ = new SymbolicSet(*X_, *U_);
        T_ = new SymbolicSet(*C_, *X2_);
        S_ = new SymbolicSet(*X_);
        Z_ = new SymbolicSet(*X_);
        solver_ = new OdeSolver(dimX, nint, tau);

        stage_ = 1;

        clog << "Initialized Compare object.\n";

    }

    /*! Destructor for a Compare object. */
    ~Compare() {
        delete X_;
        delete G_;
        delete I_;
        delete O_;
        delete X2_;
        delete U_;
        delete C_;
        delete T_;
        delete S_;
        delete Z_;
        delete solver_;
        delete Ab_;
        fclose(stderr);
    }

    /*! Saves and prints to console some information related to the reachability/always-eventually specification. */
    void saveVerifyReach() {
        cout << "X_:";
        X_->printInfo(1);
        cout << "G_:";
        G_->printInfo(1);
        cout << "I_:";
        I_->printInfo(1);
        cout << "O_:";
        O_->printInfo(1);
        cout << "U:\n";
        U_->printInfo(1);
        checkMakeDir("scots");
        X_->writeToFile("scots/X.bdd");
        G_->writeToFile("scots/G.bdd");
        I_->writeToFile("scots/I.bdd");
        O_->writeToFile("scots/O.bdd");
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

        addG(G_);
        addI(I_);
        addO(O_);

        // checking that specification is valid
        if ((G_->symbolicSet_ & O_->symbolicSet_) != ddmgr_.bddZero()) {
            error("Error: G and O have nonzero intersection.\n");
        }
        if ((I_->symbolicSet_ & O_->symbolicSet_) != ddmgr_.bddZero()) {
            error("Error: I and O have nonzero intersection.\n");
        }

        clog << "No obstacle problem with specification.\n";

        saveVerifyReach();
    }

    /*! Implementation of reachability using only a scots::FixedPoint object. For comparison with the Adaptive version.
        \param[in]  earlyBreak          If 1, the synthesis ends as soon as I meets the domain of C.
        \param[in]	verbose				If 1, prints additional information during synthesis to the log file.
    */
    void reachSCOTS(int earlyBreak, int verbose = 1) {

        TicToc tt;
        tt.tic();

        clog << "------------------------------normal SCOTS-----------------------------------\n";
        FixedPoint fp(Ab_);
        C_->symbolicSet_ = fp.reach(G_->symbolicSet_, I_->symbolicSet_, earlyBreak, verbose);
        C_->printInfo(1);
        C_->writeToFile("scots/C.bdd");
        tt.toc();
    }

    /*! Implementation of always-eventually using functions in a SCOTS::FixedPoint object. For comparison with the Adaptive version. */
    void alwaysEventuallySCOTS() {

        int curAbs = 0;
        int outerIter = 1;
        int earlyBreak = 0;
        int verbose = 1;

        TicToc tt;
        tt.tic();

        clog << "------------------------------normal SCOTS-----------------------------------\n";
        FixedPoint fp(Ab_);

        while(1) {
            clog << "Always Eventually iteration: " << outerIter << '\n';

            C_->symbolicSet_ = fp.reach(G_->symbolicSet_, I_->symbolicSet_, earlyBreak, verbose);
            Z_->symbolicSet_ = C_->symbolicSet_.ExistAbstract(U_->getCube());

            BDD preQ = fp.pre(Z_->symbolicSet_);
            BDD thing = preQ & G_->symbolicSet_;
            BDD Q = thing.ExistAbstract(X2_->getCube() & U_->getCube());

            if (G_->symbolicSet_ != Q) {
                G_->symbolicSet_ = Q;
                C_->symbolicSet_ = ddmgr_.bddZero();
                Z_->symbolicSet_ = ddmgr_.bddZero();
            }
            else {
                checkMakeDir("scots");

                SymbolicSet E(*X_);
                E.symbolicSet_ = Z_->symbolicSet_ & !Q;
                E.writeToFile("scots/E.bdd");

                SymbolicSet goalC(*C_);
                goalC.symbolicSet_ = thing;
                goalC.writeToFile("scots/goalC.bdd");

                SymbolicSet goalZ(*Z_);
                goalZ.symbolicSet_ = Q;
                goalZ.writeToFile("scots/goalZ.bdd");
                break;
            }
            outerIter += 1;
        }

        C_->printInfo(1);
        C_->writeToFile("scots/C.bdd");

        tt.toc();

    }

    /*! Saves and prints to log file some information related to the safety specification. */
    void saveVerifySafe() {
        cout << "X_:";
        X_->printInfo(1);
        cout << "S_:";
        S_->printInfo(1);
        cout << "U:\n";
        U_->printInfo(1);

        checkMakeDir("scots");
        X_->writeToFile("scots/X.bdd");
        S_->writeToFile("scots/S.bdd");
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

        addS(S_);

        saveVerifySafe();
    }

    /*! Implementation of safety using functions in a SCOTS::FixedPoint object. For comparison with the Adaptive version. */
        void safeSCOTS() {

            TicToc tt;
            tt.tic();

            clog << "------------------------------normal SCOTS-----------------------------------\n";
            FixedPoint fp(Ab_);
            C_->symbolicSet_ = fp.safe(S_->symbolicSet_, 1);
            tt.toc();

            C_->printInfo(1);
            C_->writeToFile("scots/C.bdd");
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
        Ab_ = new SymbolicModelGrowthBound<X_type, U_type>(X_, U_, X2_);
        if (readAb_ == 0) {
            Ab_->computeTransitionRelation(sysNext, radNext, X_->symbolicSet_, *solver_);
        }
        else {
            SymbolicSet T(ddmgr_, "scots/T.bdd");
            Ab_->transitionRelation_ = T.symbolicSet_;
        }

        clog << "------------------------------------------------computeAbstractions: ";
        tt.toc();

        T_->symbolicSet_ = Ab_->transitionRelation_;

        T_->printInfo(1);

        if (readAb_ == 0) {
            T_->writeToFile("scots/T.bdd");
        }

        Ab_->transitionRelation_ &= !(O_->symbolicSet_);
        clog << "Obstacles removed from transition relation.\n";

        T_->printInfo(1);

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
    void printEtaX(double* etaX, int dimX) {
        clog << "etaX: ";
        printArray(etaX, dimX);
    }

    /*! Prints information regarding the time sampling parameter to the log file. */
    void printTau(double tau) {
        clog << "tau: " << tau << '\n';

    }

};
}

#endif /* COMPARE_HH_ */
