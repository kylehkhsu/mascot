/*! \file Reach.hh
 *  Contains the Safe class. */

#ifndef SAFE_HH_
#define SAFE_HH_

#include "Adaptive.hh"

using namespace helper;

namespace scots {

/*! \class Safe
 *  \brief A class (derived from base Adaptive) that does adaptive multiscale abstraction-based synthesis for a safe specification.
 */
template<class X_type, class U_type>
class Safe: public Adaptive<X_type, U_type> {
public:
    vector<SymbolicSet*> Os_; /*!< Instance of *Xs_[i] containing unsafe (obstacle) states. */
    vector<SymbolicSet*> Ss_; /*!< Instance of *Xs_[i] containing safe states. */
    vector<SymbolicSet*> infZs_; /*!< Instance of *Xs_[i] containing projection of convergence of previous maximal fixed points. */

    /*! Constructor for a Safe object. */
    Safe(int dimX, double* lbX, double* ubX, double* etaX, double tau,
          int dimU, double* lbU, double* ubU, double* etaU,
          double* etaRatio, double tauRatio, int nint,
          int numAbs, int readXX, int readAbs, char* logFile)
        : Adaptive<X_type, U_type>(dimX, lbX, ubX, etaX, tau,
                                   dimU, lbU, ubU, etaU,
                                   etaRatio, tauRatio, nint,
                                   numAbs, readXX, readAbs, logFile)
    {
    }

    ~Safe() {
        deleteVec(Os_);
        deleteVec(Ss_);
        deleteVec(infZs_);
    }

    /*! Writes, should they exist, controllers of specified abstractions that together satisfy a safety specification. */
    void safe() {
        if (this->stage_ != 3) {
            error("Error: reach called out of order.\n");
        }

        TicToc tt;
        tt.tic();

        // removing obstacles from transition relation
        this->removeFromTs(&Os_);

        for (int curAbs = 0; curAbs < this->numAbs_; curAbs++) {
            int iter = 1;
            while (1) {
                clog << ".";
                // get pre of current abtraction's Z disjuncted with projection of converged Z from previous abstraction
                this->Cs_[curAbs]->symbolicSet_ = this->preC(this->Zs_[curAbs]->symbolicSet_ | infZs_[curAbs]->symbolicSet_, this->Ts_[curAbs]->symbolicSet_, *this->TTs_[curAbs], curAbs);
                // conjunct with safe set (maximal fixed point)
                this->Cs_[curAbs]->symbolicSet_ &= Ss_[curAbs]->symbolicSet_;
                // project onto Xs_[curAbs]
                BDD Z = this->Cs_[curAbs]->symbolicSet_.ExistAbstract(*this->notXvars_[curAbs]);

                if (Z != this->Zs_[curAbs]->symbolicSet_) { // not converged
                    this->Zs_[curAbs]->symbolicSet_ = Z; // update and continue
                }
                else { // converged
                    clog << iter << '\n';
                    this->Zs_[curAbs]->printInfo(1);
                    break;
                }
                iter += 1;
            }
            if (curAbs != this->numAbs_ - 1) {
                this->innerFinerAligned(this->Zs_[curAbs], this->infZs_[curAbs+1], curAbs); // obtain projection of converged Z onto next, finer abstraction

                this->Ts_[curAbs+1]->symbolicSet_ &= !(infZs_[curAbs+1]->symbolicSet_); // don't consider pre states that already have a controller in a coarser abstraction
                *this->TTs_[curAbs+1] &= !(infZs_[curAbs+1]->symbolicSet_); // same as above
            }
            this->Xs_[curAbs]->symbolicSet_ = this->Zs_[curAbs]->symbolicSet_ | infZs_[curAbs]->symbolicSet_; // for verification purposes
        }

        clog << "----------------------------------------safe: ";
        tt.toc();

        checkMakeDir("Z");
        saveVec(this->Zs_, "Z/Z");
        checkMakeDir("C");
        saveVec(this->Cs_, "C/C");

        cout << "Domain of all controllers:\n";
        this->Xs_[this->numAbs_-1]->printInfo(1);
        this->Xs_[this->numAbs_-1]->writeToFile("plotting/D.bdd");
    }

    /*! Initializes objects specific to the following specifications: safe.
        \param[in]	addS	Function pointer specifying the points that should be added to the potential safe set.
    */
    template<class S_type>
    void initialize(S_type addS) {
        if (this->stage_ != 1) {
            error("Error: initializeSafe called out of order.\n");
        }
        this->stage_ = 2;

        TicToc tt;
        tt.tic();
        this->initializeSpec(&Ss_, addS);
        clog << "Ss_ initialized\n";
        clog << "------------------------------------------------initializeSpec on Ss_: ";
        tt.toc();

        // empty obstacles
        for (int i = 0; i < this->numAbs_; i++) {
            SymbolicSet* O = new SymbolicSet(*this->Xs_[i]);
            Os_.push_back(O);
        }

        // maximal fixed point starts with whole set
        for (int i = 0; i < this->numAbs_; i++) {
            this->Zs_[i]->symbolicSet_ = this->ddmgr_->bddOne();
        }
        clog << "Zs_ modified to BDD-1.\n";

        for (int i = 0; i < this->numAbs_; i++) {
            SymbolicSet* infZ = new SymbolicSet(*this->Xs_[i]);
            infZs_.push_back(infZ);
        }
        clog << "infZs_ initialized to empty.\n";

        saveVerify();
    }

    /*! Saves and prints to log file some information related to the safety specification. */
    void saveVerify() {
        printVec(Ss_, "S");
        checkMakeDir("S");
        saveVec(Ss_,"S/S");
    }

    //    /*! Debugging fucntion. Implementation of basic SCOTS safety (single abstraction) in the Adaptive framework. */
    //    void safeBasicDebug(int verbose = 0) {
    //        if (numAbs_ != 1) {
    //            error("Error: For comparison with SCOTS, numAbs needs to be 1.\n");
    //        }

    //        int curAbs = 0;

    //        TicToc tt;
    //        tt.tic();

    //        SymbolicSet C(*Xs_[curAbs], *U_);

    //        int iter = 1;

    //        while (1) {
    //            // get pre(Z)
    //            C.symbolicSet_ = preC(Zs_[curAbs]->symbolicSet_, curAbs);
    //            // conjunct with safe set (maximal fixed point)
    //            C.symbolicSet_ &= Ss_[curAbs]->symbolicSet_;
    //            // project onto Xs_[curAbs]
    //            BDD Z = C.symbolicSet_.ExistAbstract(*notXvars_[curAbs]);

    //            clog << "iter: " << iter << '\n';
    //            iter += 1;

    //            if (Z != Zs_[curAbs]->symbolicSet_) {
    //                Zs_[curAbs]->symbolicSet_ = Z;
    //            }
    //            else {
    //                break;
    //            }
    //        }

    //        tt.toc();

    //        C.printInfo(1);
    //        checkMakeDir("debug");
    //        C.writeToFile("debug/C.bdd");

    //        SymbolicSet X(*Xs_[curAbs]);
    //        X.symbolicSet_ = C.symbolicSet_.ExistAbstract(U_->getCube());

    //        clog << "Domain of controller:\n";
    //        X.printInfo(1);
    //    }

};
}

#endif /* SAFE_HH_ */
